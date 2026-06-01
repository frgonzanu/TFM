import os
import re
import unicodedata
import inspect

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import silhouette_score


# ============================================================
# CONFIGURACIÓN
# ============================================================

scvi.settings.seed = 1234

input_h5ad = os.path.expanduser("~/fgonzaleznun/TFM/data/DL/sce_TPNL_scvi_in_5k.h5ad")

analysis_name = "TPNL_batch_cov_patient"
out_dir = f"results/scvi/python_outputs_{analysis_name}"
os.makedirs(out_dir, exist_ok=True)

max_epochs = 2000
n_latent = 20
gene_likelihood = "nb"

primary_label = "Tumor Primario"


# ============================================================
# FUNCIONES AUXILIARES
# ============================================================

def sanitize_filename(text):
    """
    Convierte nombres con espacios, tildes o símbolos en nombres seguros para archivos.
    """
    text = str(text)
    text = unicodedata.normalize("NFKD", text)
    text = text.encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"[^A-Za-z0-9_]+", "_", text)
    text = re.sub(r"_+", "_", text)
    text = text.strip("_")
    return text


def load_and_prepare_anndata(input_h5ad):
    """
    Carga AnnData, guarda X como layer counts, comprueba metadatos
    y convierte columnas relevantes en categóricas.
    """
    adata = sc.read_h5ad(input_h5ad)

    # scVI/LDVAE deben recibir counts crudos
    adata.layers["counts"] = adata.X.copy()

    required_cols = ["region_long", "patient"]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]

    if len(missing_cols) > 0:
        raise ValueError(f"Faltan columnas en adata.obs: {missing_cols}")

    # Convertir columnas existentes a category
    candidate_cols = ["region_long", "patient", "batch"]
    for col in candidate_cols:
        if col in adata.obs.columns:
            adata.obs[col] = adata.obs[col].astype("category")

    print("\n=============")
    print("AnnData cargado")
    print("===============")
    print(adata)

    print("\nTabla patient x region_long:")
    print(pd.crosstab(adata.obs["patient"], adata.obs["region_long"]))

    return adata



def setup_model_anndata(model_class, adata):
    """
    Registra AnnData en scvi-tools.
    """

    signature = inspect.signature(model_class.setup_anndata)
    accepts_covariates = "categorical_covariate_keys" in signature.parameters

    model_class.setup_anndata(
        adata,
        layer= "counts",
        batch_key="batch",
	categorical_covariate_keys=["patient"]
    )


def save_umaps(adata, basis, prefix, out_dir):
    """
    Guarda UMAPs
    """
    color_cols = ["region_long", "patient"]

    for col in color_cols:
        if col in adata.obs.columns:
            sc.pl.embedding(
                adata,
                basis=basis,
                color=col,
                show=False
            )
            plt.savefig(
                os.path.join(out_dir, f"{prefix}_{col}.png"),
                dpi=300,
                bbox_inches="tight"
            )
            plt.close()


def compute_silhouette_scores(adata, representation_key, output_prefix):
    """
    Calcula silhouette score para evaluar si el espacio latente separa por
    region_long o patient
    """
    X = adata.obsm[representation_key]

    rows = []

    for col in ["region_long", "patient"]:
        if col not in adata.obs.columns:
            continue

        labels = adata.obs[col].astype(str).values
        n_labels = len(np.unique(labels))

        if n_labels < 2 or n_labels >= X.shape[0]:
            continue

        try:
            score = silhouette_score(X, labels)
        except Exception as e:
            print(f"No se pudo calcular silhouette para {col}: {e}")
            score = np.nan

        rows.append({
            "representation": representation_key,
            "variable": col,
            "silhouette_score": score,
            "n_labels": n_labels
        })

    df = pd.DataFrame(rows)

    df.to_csv(
        os.path.join(out_dir, f"{output_prefix}_silhouette_scores.csv"),
        index=False
    )

    print(f"\nSilhouette scores para {representation_key}:")
    print(df)

    return df

"""
 ============================================================
 1. scVI
 ============================================================
"""

print("\n\n############################################################")
print("\nENTRENANDO scVI: TUMOR PRIMARIO VS METÁSTASIS NODO LINFÁTICO")
print("\n############################################################")

adata_scvi = load_and_prepare_anndata(input_h5ad)

setup_model_anndata(
    model_class=scvi.model.SCVI,
    adata=adata_scvi
)

model_scvi = scvi.model.SCVI(
    adata_scvi,
    n_latent=n_latent,
    gene_likelihood=gene_likelihood
)

model_scvi.train(
    max_epochs=max_epochs,
    early_stopping=True
)

adata_scvi.obsm["X_scVI"] = model_scvi.get_latent_representation()

latent_scvi = pd.DataFrame(
    adata_scvi.obsm["X_scVI"],
    index=adata_scvi.obs_names,
    columns=[f"scVI_{i+1}" for i in range(adata_scvi.obsm["X_scVI"].shape[1])]
)

latent_scvi.to_csv(
    os.path.join(out_dir, "latent_scVI_TPNL_patient_batch.csv")
)

model_scvi.save(
    os.path.join(out_dir, "model_scVI_TPNL_patient_batch"),
    overwrite=True
)

sc.pp.neighbors(
    adata_scvi,
    use_rep="X_scVI"
)

sc.tl.umap(adata_scvi)

adata_scvi.obsm["X_umap_scVI"] = adata_scvi.obsm["X_umap"].copy()

save_umaps(
    adata=adata_scvi,
    basis="umap_scVI",
    prefix="umap_scVI_TPNL_patient_batch",
    out_dir=out_dir
)

compute_silhouette_scores(
    adata=adata_scvi,
    representation_key="X_scVI",
    output_prefix="scVI_TPNL_patient_batch"
)

"""
============================================================
2. scVI: DE exploratoria por region_long
============================================================
"""

print("\n\n#########################")
print("\nDE EXPLORATORIA CON scVI")
print("\#########################")

region_categories = list(adata_scvi.obs["region_long"].cat.categories)
print("\nCategorías en region_long:")
print(region_categories)

comparisons = []

if primary_label in region_categories:
    for group in region_categories:
        if group != primary_label:
            comparisons.append((group, primary_label))
else:
    print(f"\nno se encuentra primary_label='{primary_label}' en region_long")

for group1, group2 in comparisons:
    if group1 in region_categories and group2 in region_categories:
        print(f"\nCalculando DE scVI: {group1} vs {group2}")

        de = model_scvi.differential_expression(
            adata=adata_scvi,
            groupby="region_long",
            group1=group1,
            group2=group2
        )

        safe_group1 = sanitize_filename(group1)
        safe_group2 = sanitize_filename(group2)

        de.to_csv(
            os.path.join(out_dir, f"DE_scVI_{safe_group1}_vs_{safe_group2}.csv")
        )
    else:
        print(f"\nno se encuentra: {group1} o {group2}")

"""
============================================================
3. LDVAE / LinearSCVI
============================================================
"""

print("\n\n###################################")
print("ENTRENANDO LDVAE / LinearSCVI: TPNL")
print("###################################")

adata_ldvae = load_and_prepare_anndata(input_h5ad)

setup_model_anndata(
    model_class=scvi.model.LinearSCVI,
    adata=adata_ldvae
)

model_ldvae = scvi.model.LinearSCVI(
    adata_ldvae,
    n_latent=n_latent,
    gene_likelihood=gene_likelihood
)

model_ldvae.train(
    max_epochs=max_epochs,
    early_stopping=True
)

adata_ldvae.obsm["X_LDVAE"] = model_ldvae.get_latent_representation()

latent_ldvae = pd.DataFrame(
    adata_ldvae.obsm["X_LDVAE"],
    index=adata_ldvae.obs_names,
    columns=[f"LDVAE_{i+1}" for i in range(adata_ldvae.obsm["X_LDVAE"].shape[1])]
)

latent_ldvae.to_csv(
    os.path.join(out_dir, "latent_LDVAE_TPNL_patient_batch.csv")
)

model_ldvae.save(
    os.path.join(out_dir, "model_LDVAE_TPNL_patient_batch"),
    overwrite=True
)

"""
============================================================
4. Loadings de LDVAE
============================================================
"""

print("\n\n############################################################")
print("EXTRAYENDO LOADINGS DE LDVAE")
print("############################################################")

loadings = model_ldvae.get_loadings()

if isinstance(loadings, pd.DataFrame):
    loadings_df = loadings.copy()
else:
    loadings_df = pd.DataFrame(
        loadings,
        index=adata_ldvae.var_names,
        columns=[f"LDVAE_{i+1}" for i in range(loadings.shape[1])]
    )

loadings_df.to_csv(
    os.path.join(out_dir, "LDVAE_TPNL_patient_batch_gene_loadings.csv")
)

top_genes_all = []

for factor in loadings_df.columns:
    tmp = loadings_df[[factor]].copy()
    tmp["gene"] = tmp.index
    tmp["abs_loading"] = tmp[factor].abs()
    tmp = tmp.sort_values("abs_loading", ascending=False)
    tmp["factor"] = factor
    top_genes_all.append(tmp.head(50))

top_genes_ldvae = pd.concat(top_genes_all, axis=0)

top_genes_ldvae.to_csv(
    os.path.join(out_dir, "LDVAE_TPNL_patient_batch_top50_genes_per_factor.csv"),
    index=False
)

"""
============================================================
5. Top genes positivos y negativos por factor LDVAE
============================================================
"""

for factor in loadings_df.columns:
    tmp = loadings_df[[factor]].copy()
    tmp["gene"] = tmp.index
    tmp["abs_loading"] = tmp[factor].abs()

    top_pos = tmp.sort_values(factor, ascending=False).head(50)
    top_neg = tmp.sort_values(factor, ascending=True).head(50)
    top_abs = tmp.sort_values("abs_loading", ascending=False).head(100)

    safe_factor = sanitize_filename(factor)

    top_pos.to_csv(
        os.path.join(out_dir, f"{safe_factor}_top50_positive_loadings.csv"),
        index=False
    )

    top_neg.to_csv(
        os.path.join(out_dir, f"{safe_factor}_top50_negative_loadings.csv"),
        index=False
    )

    top_abs.to_csv(
        os.path.join(out_dir, f"{safe_factor}_top100_abs_loadings.csv"),
        index=False
    )

"""
============================================================
6. Relación factores LDVAE - tejido
============================================================
"""

print("\n\n################################")
print("RELACIÓN FACTORES LDVAE - TEJIDO")
print("################################")

ldvae_df = pd.DataFrame(
    adata_ldvae.obsm["X_LDVAE"],
    index=adata_ldvae.obs_names,
    columns=[f"LDVAE_{i+1}" for i in range(adata_ldvae.obsm["X_LDVAE"].shape[1])]
)

metadata_cols = ["region_long", "patient"]

for col in metadata_cols:
    ldvae_df[col] = adata_ldvae.obs[col].values

ldvae_df.to_csv(
    os.path.join(out_dir, "LDVAE_TPNL_patient_batch_factors_with_metadata.csv")
)

factor_cols = [c for c in ldvae_df.columns if c.startswith("LDVAE_")]

summary_by_region = ldvae_df.groupby("region_long")[factor_cols].agg(
    ["mean", "median", "std"]
)

summary_by_region.to_csv(
    os.path.join(out_dir, "LDVAE_TPNL_patient_batch_factor_summary_by_region.csv")
)

summary_by_patient_region = (
    ldvae_df
    .groupby(["patient", "region_long"])[factor_cols]
    .mean()
    .reset_index()
)

summary_by_patient_region.to_csv(
    os.path.join(out_dir, "LDVAE_TPNL_patient_batch_factor_means_by_patient_region.csv"),
    index=False
)

"""
============================================================
7. Boxplots LDVAE por región
============================================================
"""

region_order = [
    "Tumor Primario",
    "Metastasis Nodo Linfático"
]

region_order = [x for x in region_order if x in ldvae_df["region_long"].unique()]

for factor in factor_cols:
    plt.figure(figsize=(6, 4))

    sns.boxplot(
        data=ldvae_df,
        x="region_long",
        y=factor,
        order=region_order if len(region_order) > 0 else None
    )

    sns.stripplot(
        data=ldvae_df,
        x="region_long",
        y=factor,
        order=region_order if len(region_order) > 0 else None,
        color="black",
        alpha=0.4,
        size=2
    )

    plt.xticks(rotation=45, ha="right")
    plt.title(f"{factor} by tissue region")
    plt.tight_layout()

    safe_factor = sanitize_filename(factor)

    plt.savefig(
        os.path.join(out_dir, f"boxplot_TPNL_patient_batch_{safe_factor}_by_region.png"),
        dpi=300
    )

    plt.close()

"""
============================================================
8. UMAP LDVAE
============================================================
"""

print("\n\n##########")
print("UMAP LDVAE")
print("##########")

sc.pp.neighbors(
    adata_ldvae,
    use_rep="X_LDVAE",
    key_added="neighbors_LDVAE"
)

sc.tl.umap(
    adata_ldvae,
    neighbors_key="neighbors_LDVAE"
)

adata_ldvae.obsm["X_umap_LDVAE"] = adata_ldvae.obsm["X_umap"].copy()

save_umaps(
    adata=adata_ldvae,
    basis="umap_LDVAE",
    prefix="umap_LDVAE_TPNL_patient_batch",
    out_dir=out_dir
)

compute_silhouette_scores(
    adata=adata_ldvae,
    representation_key="X_LDVAE",
    output_prefix="LDVAE_TPNL_patient_batch"
)

"""
============================================================
9. Guardar AnnData finales
============================================================
"""

adata_scvi.write_h5ad(
    os.path.join(out_dir, "adata_scVI_TPNL_patient_bact_results.h5ad")
)

adata_ldvae.write_h5ad(
    os.path.join(out_dir, "adata_LDVAE_TPNL_patient_batch_results.h5ad")
)

print("\n\n####")
print("DONE")
print("####")
