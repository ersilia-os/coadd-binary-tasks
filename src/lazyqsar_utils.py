import json
import os

import numpy as np
from lazyqsar.qsar import LazyClassifierQSAR
from lazyqsar.utils.logging import logger
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit


def run_crossval(smiles_list, y, n_folds, output_path,
                 mode="slow", test_size=0.20, random_state=42):
    """Stratified shuffle CV saving predict_proba, predict_score and predict_rank per fold.

    Parameters
    ----------
    smiles_list  : list[str]
    y            : array-like of int (binary labels)
    n_folds      : int
    output_path  : str  — path to write the JSON report
    mode         : str  — LazyBinaryQSAR mode ("fast" or "slow"), default "fast"
    test_size    : float (default 0.20)
    random_state : int   (default 42)

    Returns
    -------
    dict  — report keyed by fold index (str)
    """
    smiles_arr = np.array(smiles_list)
    y_arr = np.array(y)
    sss = StratifiedShuffleSplit(
        n_splits=n_folds, test_size=test_size, random_state=random_state
    )
    report = {}
    for fold, (train_idx, test_idx) in enumerate(sss.split(smiles_arr, y_arr)):
        smiles_train = smiles_arr[train_idx].tolist()
        smiles_test = smiles_arr[test_idx].tolist()
        y_train = y_arr[train_idx]
        y_test = y_arr[test_idx]

        model = LazyClassifierQSAR(mode=mode)
        model.fit(smiles_list=smiles_train, y=y_train)

        y_hat   = model.predict_proba(smiles_list=smiles_test)[:, 1]
        y_score = model.predict_score(smiles_list=smiles_test)[:, 1]
        y_rank  = model.predict_rank(smiles_list=smiles_test)[:, 1]

        roc_auc = float(roc_auc_score(y_test, y_hat))
        logger.info(f"Fold {fold} ROC-AUC: {roc_auc:.4f}")

        report[str(fold)] = {
            "y_true":  y_test.tolist(),
            "y_hat":   y_hat.tolist(),
            "y_score": y_score.tolist(),
            "y_rank":  y_rank.tolist(),
            "roc_auc": roc_auc,
        }

    with open(output_path, "w") as f:
        json.dump(report, f, indent=2)
    return report

def predict_and_save(model, smiles_list, output_path, y_true=None):
    """Call predict_proba, predict_score and predict_rank; save all to JSON.

    Parameters
    ----------
    model       : fitted LazyClassifierQSAR
    smiles_list : list[str]
    output_path : str
    y_true      : array-like of int, optional

    Returns
    -------
    np.ndarray — predict_proba[:, 1]
    """
    y_hat   = model.predict_proba(smiles_list=smiles_list)[:, 1]
    y_score = model.predict_score(smiles_list=smiles_list)[:, 1]
    y_rank  = model.predict_rank(smiles_list=smiles_list)[:, 1]

    result = {
        "smiles":  smiles_list,
        "y_hat":   y_hat.tolist(),
        "y_score": y_score.tolist(),
        "y_rank":  y_rank.tolist(),
    }
    if y_true is not None:
        result["y_true"] = [int(v) for v in y_true]

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)

    return y_hat