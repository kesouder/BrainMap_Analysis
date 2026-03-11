import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
def subgroup_connectivity(corr_df, labels, label_index=None):
    """
    corr_df: (n x n) DataFrame
    labels: array-like length n (in same order as corr_df.index unless label_index provided)
    label_index: optional index aligned to corr_df.index
    """
    if label_index is None:
        label_index = corr_df.index

    lab = pd.Series(labels, index=label_index).reindex(corr_df.index)
    C = corr_df.values
    n = C.shape[0]

    within = []
    between = []
    for i in range(n):
        for j in range(i+1, n):
            if lab.iloc[i] == lab.iloc[j]:
                within.append(C[i,j])
            else:
                between.append(C[i,j])

    return {
        "within_mean": float(np.mean(within)) if within else np.nan,
        "between_mean": float(np.mean(between)) if between else np.nan,
        "within_vals": np.array(within),
        "between_vals": np.array(between),
        "labels": lab
    }

def plot_group_corr(res):
    plt.figure(figsize=(6,4))
    plt.hist(res["within_vals"], alpha=0.6, label="within")
    plt.hist(res["between_vals"], alpha=0.6, label="between", color="r")
    plt.legend()
    plt.xlabel("correlation")
    plt.ylabel("count")
    plt.savefig("images/groupwise_hist.jpg")
    plt.show()