import matplotlib.patches as mpatches
import numpy as np
import stylia as st
from stylia import DivergingColormap, FadingColormap, SpectralColormap
from sklearn.metrics import roc_curve, auc
import json
import pandas as pd

st.set_format("print")
st.set_style("ersilia")
nc = st.NamedColors()


def load_report(model_path, dataset, cutoff):
    """Load cross-validation report from JSON. Returns the raw dict keyed by fold.
    When cutoff is None, loads report_{dataset}.json (no cutoff suffix)."""
    if cutoff is None:
        path = f"{model_path}/report_{dataset}.json"
    else:
        path = f"{model_path}/report_{dataset}_{cutoff}um.json"
    with open(path) as f:
        return json.load(f)


def _folds_from_report(data):
    """Convert report dict to list of (y_true, y_hat) arrays."""
    return [
        (np.asarray(v["y_true"], dtype=float), np.asarray(v["y_hat"], dtype=float))
        for v in data.values()
    ]


def plot_roc_mean(ax, data, dataset, cutoff):
    """
    Plot mean ROC curve ± 1 SD band across cross-validation folds.

    Parameters
    ----------
    ax : matplotlib Axes
    data : dict  — raw report dict (keyed by fold)
    dataset : str
    cutoff : int  — µM cutoff used for binarization
    """
    folds = _folds_from_report(data)
    mean_fpr = np.linspace(0, 1, 200)
    tprs, aucs = [], []

    for y_true, y_hat in folds:
        mask = np.isfinite(y_hat) & np.isfinite(y_true)
        y_true, y_hat = y_true[mask], y_hat[mask]
        fpr, tpr, _ = roc_curve(y_true, y_hat)
        tpr_interp = np.interp(mean_fpr, fpr, tpr, left=0.0, right=1.0)
        tpr_interp[0] = 0.0
        tprs.append(tpr_interp)
        aucs.append(auc(fpr, tpr))

    tprs = np.vstack(tprs)
    mean_tpr = tprs.mean(axis=0)
    std_tpr = tprs.std(axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)

    ax.plot(mean_fpr, mean_tpr, color=nc.plum,
            label=f"Mean ROC (AUC = {mean_auc:.2f} ± {std_auc:.2f})")
    ax.fill_between(
        mean_fpr,
        np.maximum(mean_tpr - std_tpr, 0),
        np.minimum(mean_tpr + std_tpr, 1),
        color=nc.purple, alpha=0.25, label="±1 SD",
    )
    ax.plot([0, 1], [0, 1], "--", color=nc.gray)
    st.label(ax,
             xlabel="False Positive Rate",
             ylabel="True Positive Rate",
             title=f"ROC — {dataset}" + (f" {cutoff}µM" if cutoff is not None else ""))
    ax.legend()


def plot_swarm(ax, data, dataset, cutoff):
    """
    Jittered scatter of predicted scores split by true label (first fold only).

    Parameters
    ----------
    ax : matplotlib Axes
    data : dict  — raw report dict (keyed by fold)
    dataset : str
    cutoff : int  — µM cutoff used for binarization
    """
    folds = _folds_from_report(data)
    y_true, y_hat = folds[0]
    y_true = y_true.astype(int)

    y0 = y_hat[y_true == 0]
    y1 = y_hat[y_true == 1]
    rng = np.random.default_rng(42)
    rng.shuffle(y0)
    rng.shuffle(y1)

    def jitter(n, center):
        return center + rng.uniform(-0.12, 0.12, size=n)

    ax.scatter(jitter(len(y0), 0), y0, alpha=0.5, s=st.get_markersize("small"),
               color=nc.gray, edgecolors="none", label="Neg")
    ax.scatter(jitter(len(y1), 1), y1, alpha=0.5, s=st.get_markersize("small"),
               color=nc.orange, edgecolors="none", label="Pos")
    ax.set_xlim(-0.5, 1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Neg", "Pos"])
    ax.set_ylim(-0.05, 1.05)
    st.label(ax,
             xlabel="Test labels",
             ylabel="Score",
             title=f"Scores — {dataset}" + (f" {cutoff}µM" if cutoff is not None else ""))
    ax.legend()


def plot_scores(ax, data, title, fold=0):
    """
    4-box grouped boxplot: probability and rank for neg/pos from one fold.

    Two groups: [proba neg, proba pos] and [rank neg, rank pos].
    Neg boxes are gray, pos boxes are orange.
    Box lines, whiskers, caps and median are drawn in plum.

    Parameters
    ----------
    ax    : matplotlib Axes
    data  : dict  — raw report dict (keyed by fold index as str)
    title : str
    fold  : int   — fold index to display (default 0)
    """
    v = data[str(fold)]
    y_true = np.asarray(v["y_true"], dtype=int)
    y_hat  = np.asarray(v["y_hat"],  dtype=float)
    y_rank = np.asarray(v["y_rank"], dtype=float)

    proba_neg = y_hat[y_true == 0]
    proba_pos = y_hat[y_true == 1]
    rank_neg  = y_rank[y_true == 0]
    rank_pos  = y_rank[y_true == 1]

    positions = [0, 1, 3, 4]
    bplot = ax.boxplot(
        [proba_neg, proba_pos, rank_neg, rank_pos],
        positions=positions,
        patch_artist=True,
        medianprops=dict(color=nc.plum, linewidth=0.5),
        boxprops=dict(color=nc.plum, linewidth=0.5),
        whiskerprops=dict(color=nc.plum, linewidth=0.5),
        capprops=dict(color=nc.plum, linewidth=0.5),
        flierprops=dict(marker="o", markeredgecolor=nc.plum,
                        markerfacecolor="none", markersize=2, markeredgewidth=0.5),
        manage_ticks=False,
    )
    colors = [nc.gray, nc.orange, nc.gray, nc.orange]
    for patch, color in zip(bplot["boxes"], colors):
        patch.set_facecolor(color)

    ax.set_xticks([0, 1, 3, 4])
    ax.set_xticklabels(["Neg", "Pos", "Neg", "Pos"])
    ax.set_xlim(-0.8, 5.3)
    st.label(ax, ylabel="proba / rank", xlabel="Proba vs Rank", title=title)


def plot_class_balance(ax, data, title=""):
    """Bar chart of Pos/Neg class counts from all CV folds combined.

    Parameters
    ----------
    ax    : matplotlib Axes
    data  : dict  — raw report dict (keyed by fold)
    title : str
    """
    all_y_true = []
    for v in data.values():
        all_y_true.extend(v["y_true"])
    all_y_true = np.array(all_y_true, dtype=int)
    n_pos = int(all_y_true.sum())
    n_neg = int((all_y_true == 0).sum())

    bars = ax.bar([0, 1], [n_neg, n_pos], color=[nc.gray, nc.orange])
    for bar, count in zip(bars, [n_neg, n_pos]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                str(count), ha="center", va="bottom")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Neg", "Pos"])
    ax.set_ylim(0, max(n_neg, n_pos) * 1.1)
    st.label(ax, ylabel="Count", xlabel="", title=title)


def plot_roc_folds(ax, data):
    """
    Plot each fold's ROC curve colored by its AUROC.
    Annotates mean AUROC and class counts.

    Parameters
    ----------
    ax : matplotlib Axes
    data : dict  — raw report dict (keyed by fold)
    """
    cmap = FadingColormap("plum")
    auroc_values = [v["roc_auc"] for v in data.values()]
    cmap.fit([0.5, 1])

    for v in data.values():
        color = cmap.transform([v["roc_auc"]])[0]
        fpr, tpr, _ = roc_curve(
            np.array(v["y_true"]), np.array(v["y_hat"])
        )
        ax.plot(fpr, tpr, color=color)

    mean_auroc = np.mean(auroc_values)
    std_auc    = np.std(auroc_values)

    ax.plot([], [], linestyle="none",
            label=f"Mean AUC={mean_auroc:.3f}±{std_auc:.2f}")
    ax.legend(loc="lower right")


def plot_corr_matrix(ax, corr, title, show_labels=True):
    """Annotated Spearman correlation heatmap using DivergingColormap plum_mint.

    Parameters
    ----------
    ax          : matplotlib Axes
    corr        : pd.DataFrame — square correlation matrix
    title       : str
    show_labels : bool — show tick labels and cell annotations (default True)
    """
    cmap = DivergingColormap("plum_mint")
    cmap.fit([-1, 1])
    n = len(corr)
    labels = corr.columns.tolist()

    for i in range(n):
        for j in range(n):
            val = corr.iloc[i, j]
            color = cmap.transform([val])[0]
            ax.add_patch(mpatches.Rectangle(
                (j, n - i - 1), 1, 1, color=color, ec="white"
            ))
            if show_labels:
                ax.text(j + 0.5, n - i - 0.5, f"{val:.2f}",
                        ha="center", va="center",
                        color="white" if abs(val) > 0.5 else nc.plum)

    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    if show_labels:
        ax.set_xticks([i + 0.5 for i in range(n)])
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticks([i + 0.5 for i in range(n)])
        ax.set_yticklabels(list(reversed(labels)))
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    st.label(ax, xlabel="", ylabel="", title=title)

def compute_hit_overlap(pred_df, k):
    """Pairwise hit-overlap matrix: overlap[i,j] = |top_k(i) ∩ top_k(j)| / k."""
    k = min(k, len(pred_df))
    models = pred_df.columns.tolist()
    top_k = {m: set(pred_df[m].nlargest(k).index) for m in models}
    n = len(models)
    mat = np.zeros((n, n))
    for i, mi in enumerate(models):
        for j, mj in enumerate(models):
            mat[i, j] = len(top_k[mi] & top_k[mj]) / k
    return pd.DataFrame(mat, index=models, columns=models)

def plot_overlap_matrix(ax, overlap_df, title, show_labels=True):
    """Annotated hit-overlap heatmap using DivergingColormap (range 0–1)."""
    cmap = DivergingColormap("purple_orange")
    cmap.fit([0, 1])
    n = len(overlap_df)
    labels = overlap_df.columns.tolist()

    for i in range(n):
        for j in range(n):
            val = float(overlap_df.iloc[i, j])
            color = cmap.transform([val])[0]
            ax.add_patch(mpatches.Rectangle(
                (j, n - i - 1), 1, 1, color=color, ec="white"
            ))
            if show_labels:
                ax.text(j + 0.5, n - i - 0.5, f"{val:.2f}",
                        ha="center", va="center",
                        color=nc.plum)

    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    if show_labels:
        ax.set_xticks([i + 0.5 for i in range(n)])
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticks([i + 0.5 for i in range(n)])
        ax.set_yticklabels(list(reversed(labels)))
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    st.label(ax, xlabel="", ylabel="", title=title)


def plot_cross_corr(ax, corr_df, title):
    """Heatmap for a rectangular cross-correlation matrix (no cell annotations).

    corr_df rows = one model group (e.g. citotox), cols = another (e.g. bioactivity).
    Plotted with rows on y-axis, cols on x-axis.
    """
    cmap = DivergingColormap("plum_mint")
    cmap.fit([-1, 1])
    n_rows, n_cols = corr_df.shape
    row_labels = corr_df.index.tolist()
    col_labels  = corr_df.columns.tolist()

    for i in range(n_rows):
        for j in range(n_cols):
            val = float(corr_df.iloc[i, j])
            color = cmap.transform([val])[0]
            ax.add_patch(mpatches.Rectangle(
                (j, n_rows - i - 1), 1, 1, color=color, ec="white"
            ))

    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.set_xticks([j + 0.5 for j in range(n_cols)])
    ax.set_xticklabels(col_labels, rotation=90, ha="center", fontsize=5)
    ax.set_yticks([i + 0.5 for i in range(n_rows)])
    ax.set_yticklabels(list(reversed(row_labels)), fontsize=7)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    st.label(ax, xlabel="", ylabel="", title=title)
