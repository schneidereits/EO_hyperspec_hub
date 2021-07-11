# =====================================================================================================================
# Linear spectral mixing to create synthetic training data
# =====================================================================================================================
import numpy as np
import pandas as pd


# linear-spectral-mixing analysis function
def lsma(df, col_response=None, response_id=None, n=1000, within_class_mixture=True, response_mixture=False,
         includeEndmember=True, targetRange=(0, 1), mix_complexity=None, p_mix_complexity=None):
    """
    Linear Spectral Mixture Analysis function to create synthetic training data from endmembers.

    :param df: Dataframe containing input features to be mixed as well as the response variable.
    :param col_response: (string) Column name of response variable.
    :param response_id: (int) Numeric value corresponding to target class of "col_response".
    :param n: (int) Number of synthetic features to create.
    :param within_class_mixture: (bool) Allow mixtures within classes apart from target class.
    :param response_mixture: (bool) Allow mixtures within the target class.
    :param includeEndmember: (bool) Include input endmembers in output.
    :param targetRange: (int, tuple) Tuple of boundary values of the desired target range.
    :param mix_complexity: (int, list) List of integers referring to number of classes to be mixed. E.g. [2, 3] means
                                       that there will be mixtures of 2 and 3 classes
    :param p_mix_complexity: (float, list) List of floats referring to the probabilities associated with the
                                           mix_complexity, hence the expected frequency of certain mixtures.
    :return: Dataframe with synthetic mixtures of predictor and response variable.
    """

    if mix_complexity is None:
        mix_complexity = [2, 3, 4]
    if p_mix_complexity is None:
        p_mix_complexity = [0.7, 0.2, 0.1]

    response = np.asarray(df[col_response])
    unique_response = np.unique(response)
    classes = len(unique_response)
    features = np.asarray(df.drop([col_response], axis=1)).T  # bands in rows features in columns

    classLikelihoods = {i + 1: len(np.where(response == i + 1)[0]) / len(response) for i in range(classes)}

    # cache label indices and setup 0%/100% fractions from class labels
    indices = dict()
    zeroOneFractions = np.zeros((classes, features.shape[1]), dtype=np.float32)
    for label in range(1, classes + 1):
        indices[label] = np.where(response == label)[0]
        zeroOneFractions[label - 1, indices[label]] = 1.

    # create mixtures
    mixtures = list()
    fractions = list()

    classLikelihoods2 = {k: v / (1 - classLikelihoods[response_id]) for k, v in classLikelihoods.items() if k != response_id}

    for i in range(n):

        # get mixing complexity
        complexity = np.random.choice(mix_complexity, p=p_mix_complexity)

        # define current target class
        l_response = [response_id]

        # ...
        if within_class_mixture:
            if response_mixture:
                l_response.extend(np.random.choice(list(classLikelihoods.keys()), size=complexity - 1, replace=True,
                                                p=list(classLikelihoods.values())))
            else:
                l_response.extend(np.random.choice(list(classLikelihoods2.keys()), size=complexity - 1, replace=True,
                                                   p=list(classLikelihoods2.values())))
        else:
            l_response.extend(np.random.choice(list(classLikelihoods2.keys()), size=complexity - 1, replace=False,
                                                p=list(classLikelihoods2.values())))

        drawnIndices = [np.random.choice(indices[label]) for label in l_response]
        drawnFeatures = features[:, drawnIndices]
        drawnFractions = zeroOneFractions[:, drawnIndices]

        randomWeights = list()
        for i in range(complexity - 1):
            if i == 0:
                weight = np.random.random() * (targetRange[1] - targetRange[0]) + targetRange[0]
            else:
                weight = np.random.random() * (1. - sum(randomWeights))
            randomWeights.append(weight)
        randomWeights.append(1. - sum(randomWeights))

        assert sum(randomWeights) == 1.
        mixtures.append(np.sum(drawnFeatures * randomWeights, axis=1))
        fractions.append(np.sum(drawnFractions * randomWeights, axis=1)[response_id - 1])

    if includeEndmember:
        mixtures.extend(features.T)
        fractions.extend(np.float32(response == response_id))  # 1. for target class, 0. for the rest

    # convert to df
    df_final = pd.DataFrame(np.column_stack([np.repeat(response_id, len(mixtures)), mixtures, fractions]),
                            columns=list(df.columns)+['fraction'])
    return df_final


# input
input_csv = ""
output_csv = ""  # string with .csv ending; file does not need to exist
df = pd.read_csv(input_csv)  # .csv table
'''
- columns: one column holding class_id as integer (e.g., 1, 2, ..., n), the remaining columns are bands
- each row represents a single pure endmember point
- cleaned of nodata values, only valid observations (otherwise they might be mixed in)
'''
target_attr = ''  # name of column which holds the class_id
n_samples = 2500  # number of synthetically mixed training points to be generated

# run
unique_classes = np.unique(df[target_attr])  # retrieved the unique classes to mix n_samples for each target class
df_fraction = pd.DataFrame()
for i in unique_classes:
    df_fraction = df_fraction.append(lsma(df, col_response=target_attr,
                                          response_id=i, n=n_samples, mix_complexity=[2, 3, 4],
                                          p_mix_complexity=[0.75, 0.20, 0.05], targetRange=(0, 1),
                                          within_class_mixture=True, response_mixture=True, includeEndmember=True))

df_fraction[target_attr] = df_fraction[target_attr].astype('int')
df_fraction.to_csv(output_csv, index=False)


# EOF
