import pandas as pd
import hashlib
import json

B = pd.read_csv(r"D:\rna-adt-transportability-boundary\supplementary\B_fold_level_model_table.csv")

# Compute a reproducible per-fold SHA-256 over the full config + observed/predicted.
# This makes each fold's record independently verifiable (Issue 5: replay fold-by-fold).
def fold_hash(row):
    payload = {
        "held_out_target": row["held_out_target"],
        "training_targets": row["training_targets"],
        "selected_features": row["selected_features"],
        "scaler_mean": row["scaler_mean"],
        "scaler_scale": row["scaler_scale"],
        "ridge_alpha": row["ridge_alpha"],
        "ridge_solver": row["ridge_solver"],
        "observed_value": row["observed_value"],
        "predicted_value": row["predicted_value"],
    }
    s = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(s.encode("utf-8")).hexdigest()

B["fold_sha256"] = B.apply(fold_hash, axis=1)

# reorder: put hash right after fold_id
cols = list(B.columns)
cols.remove("fold_sha256")
cols.insert(1, "fold_sha256")
B = B[cols]

B.to_csv(r"D:\rna-adt-transportability-boundary\supplementary\B_fold_level_model_table.csv", index=False)

print(B[["fold_id","fold_sha256","held_out_target","ridge_alpha","observed_value","predicted_value"]].to_string(index=False))
print()
print("rows:", len(B), "| unique hashes:", B["fold_sha256"].nunique())
