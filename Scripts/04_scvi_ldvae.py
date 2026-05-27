import os
import re
import unicodedata

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

base_out_dir = "results/scvi/TPNL_no_batch_two_explorations"
os.makedirs(base_out_dir, exist_ok=True)

max_epochs = 2000
n_latent = 20
gene_likelihood = "nb"

primary_label = "Tumor Primario"

explorations = {
    "no_correction": {
        "description": "Sin corrección explícita por paciente",
        "batch_key": None
    },
    "patient_corrected": {
        "description": "Paciente usado como batch_key",
        "batch_key": "patient"
    }
}


# ============================================================
# FUNCIONES AUXILIARES
# ============================================================

def sanitize_filename(text):
    text = str(text)
    text = unicodedata.normalize("NFKD", text)
    text = text.encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"[^A-Za-z0-9_]+", "_", text)
    text = re.sub(r"_+", "_", text)
    text = text.strip("_")
    return text


def load_and_prepare_anndata(input_h5ad):
    """
    Carga AnnData, guarda X como layer counts y comprueba metadatos.
    En este análisis no se usa batch porque CRC01/CRC02 han sido excluidos.
    """
    adata = sc.read_h5ad(input_h5ad)

    adata.layers["counts"] = adata.X.copy()

    required_cols = ["region_long", "patient"]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]

    if len(missing_cols) > 0:
        raise ValueError(f"Faltan columnas en adata.obs: {missing_cols}")

    for col in ["region_long", "patient"]:
        adata.obs[col] = adata.obs[col].astype("category")

    print("\n==============================")
    print("AnnData cargado")
    print("==============================")
    print(adata)

    print("\nDistribución region_long:")
    print(adata.obs["region_long"].value_counts())

    print("\nDistribución patient:")
    print(adata.obs["patient"].value_counts())

    print("\nTabla patient x region_long:")
    print(pd.crosstab(adata.obs["patient"], adata.obs["region_long"]))

    return adata


def setup_model_anndata(model_class, adata, batch_key=None):
    """
    Registra AnnData en scvi-tools.

    batch_key=None:
        modelo sin corrección explícita.

    batch_key="patient":
        paciente usado como factor de ajuste.
    """
    setup_kwargs = {
        "layer": "counts"
    }

    if batch_key is not None:
        if batch_key not in adata.obs.columns:
            raise ValueError(f"No existe la columna '{batch_key}' en adata.obs.")

        adata.obs[batch_key] = adata.obs[batch_key].astype("category")
        setup_kwargs["batch_key"] = batch_key

    print(f"\nRegistrando {model_class.__name__} con:")
    print(setup_kwargs)

    model_class.setup_anndata(
        adata,
        **setup_kwargs
    )


def save_umaps(adata, basis, prefix, out_dir):
    """
    Guarda UMAPs coloreados por región y paciente.
    """
    for col in ["region_long", "patient"]:
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


def compute_silhouette_scores(adata, representation_key, output_prefix, out_dir):
    """
    Calcula silhouette score para ver cuánto separa el espacio latente
    por región o por paciente.
    """
    X = adata.obsm[representation_key]

    rows = []

    for col in ["region_long", "patient"]:
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


# ============================================================
# scVI
# ============================================================

def run_scvi(input_h5ad, exploration_name, exploration_config, out_dir):
    print("\n\n############################################################")
    print(f"ENTRENANDO scVI: {exploration_name}")
    print(exploration_config["description"])
    print("############################################################")

    adata = load_and_prepare_anndata(input_h5ad)

    setup_model_anndata(
        model_class=scvi.model.SCVI,
        adata=adata,
        batch_key=exploration_config["batch_key"]
    )

    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        gene_likelihood=gene_likelihood
    )

    model.train(
        max_epochs=max_epochs,
        early_stopping=True
    )

    adata.obsm["X_scVI"] = model.get_latent_representation()

    latent = pd.DataFrame(
        adata.obsm["X_scVI"],
        index=adata.obs_names,
        columns=[f"scVI_{i+1}" for i in range(adata.obsm["X_scVI"].shape[1])]
    )

    latent.to_csv(
        os.path.join(out_dir, f"latent_scVI_{exploration_name}.csv")
    )

    model.save(
        os.path.join(out_dir, f"model_scVI_{exploration_name}"),
        overwrite=True
    )

    sc.pp.neighbors(
        adata,
        use_rep="X_scVI"
    )

    sc.tl.umap(adata)

    adata.obsm["X_umap_scVI"] = adata.obsm["X_umap"].copy()

    save_umaps(
        adata=adata,
        basis="umap_scVI",
        prefix=f"umap_scVI_{exploration_name}",
        out_dir=out_dir
    )

    compute_silhouette_scores(
        adata=adata,
        representation_key="X_scVI",
        output_prefix=f"scVI_{exploration_name}",
        out_dir=out_dir
    )

    run_scvi_de(
        adata=adata,
        model=model,
        exploration_name=exploration_name,
        out_dir=out_dir
    )

    adata.write_h5ad(
        os.path.join(out_dir, f"adata_scVI_{exploration_name}_results.h5ad")
    )

    return adata, model


def run_scvi_de(adata, model, exploration_name, out_dir):
    print("\n\n#########################")
    print(f"DE EXPLORATORIA CON scVI: {exploration_name}")
    print("#########################")

    region_categories = list(adata.obs["region_long"].cat.categories)

    print("\nCategorías en region_long:")
    print(region_categories)

    if primary_label not in region_categories:
        print(f"\nNo se encuentra primary_label='{primary_label}' en region_long")
        return

    comparisons = [
        (group, primary_label)
        for group in region_categories
        if group != primary_label
    ]

    for group1, group2 in comparisons:
        print(f"\nCalculando DE scVI: {group1} vs {group2}")

        de = model.differential_expression(
            adata=adata,
            groupby="region_long",
            group1=group1,
            group2=group2
        )

        safe_group1 = sanitize_filename(group1)
        safe_group2 = sanitize_filename(group2)

        de.to_csv(
            os.path.join(
                out_dir,
                f"DE_scVI_{exploration_name}_{safe_group1}_vs_{safe_group2}.csv"
            )
        )


# ============================================================
# LDVAE / LinearSCVI
# ============================================================

def run_ldvae(input_h5ad, exploration_name, exploration_config, out_dir):
    print("\n\n############################################################")
    print(f"ENTRENANDO LDVAE / LinearSCVI: {exploration_name}")
    print(exploration_config["description"])
    print("############################################################")

    adata = load_and_prepare_anndata(input_h5ad)

    setup_model_anndata(
        model_class=scvi.model.LinearSCVI,
        adata=adata,
        batch_key=exploration_config["batch_key"]
    )

    model = scvi.model.LinearSCVI(
        adata,
        n_latent=n_latent,
        gene_likelihood=gene_likelihood
    )

    model.train(
        max_epochs=max_epochs,
        early_stopping=True
    )

    adata.obsm["X_LDVAE"] = model.get_latent_representation()

    latent = pd.DataFrame(
        adata.obsm["X_LDVAE"],
        index=adata.obs_names,
        columns=[f"LDVAE_{i+1}" for i in range(adata.obsm["X_LDVAE"].shape[1])]
    )

    latent.to_csv(
        os.path.join(out_dir, f"latent_LDVAE_{exploration_name}.csv")
    )

    model.save(
        os.path.join(out_dir, f"model_LDVAE_{exploration_name}"),
        overwrite=True
    )

    extract_ldvae_loadings(
        adata=adata,
        model=model,
        exploration_name=exploration_name,
        out_dir=out_dir
    )

    summarize_ldvae_factors(
        adata=adata,
        exploration_name=exploration_name,
        out_dir=out_dir
    )

    sc.pp.neighbors(
        adata,
        use_rep="X_LDVAE",
        key_added="neighbors_LDVAE"
    )

    sc.tl.umap(
        adata,
        neighbors_key="neighbors_LDVAE"
    )

    adata.obsm["X_umap_LDVAE"] = adata.obsm["X_umap"].copy()

    save_umaps(
        adata=adata,
        basis="umap_LDVAE",
        prefix=f"umap_LDVAE_{exploration_name}",
        out_dir=out_dir
    )

    compute_silhouette_scores(
        adata=adata,
        representation_key="X_LDVAE",
        output_prefix=f"LDVAE_{exploration_name}",
        out_dir=out_dir
    )

    adata.write_h5ad(
        os.path.join(out_dir, f"adata_LDVAE_{exploration_name}_results.h5ad")
    )

    return adata, model


def extract_ldvae_loadings(adata, model, exploration_name, out_dir):
    print("\n\n############################################################")
    print(f"EXTRAYENDO LOADINGS DE LDVAE: {exploration_name}")
    print("############################################################")

    loadings = model.get_loadings()

    if isinstance(loadings, pd.DataFrame):
        loadings_df = loadings.copy()
    else:
        loadings_df = pd.DataFrame(
            loadings,
            index=adata.var_names,
            columns=[f"LDVAE_{i+1}" for i in range(loadings.shape[1])]
        )

    loadings_df.to_csv(
        os.path.join(out_dir, f"LDVAE_{exploration_name}_gene_loadings.csv")
    )

    top_genes_all = []

    for factor in loadings_df.columns:
        tmp = loadings_df[[factor]].copy()
        tmp["gene"] = tmp.index
        tmp["abs_loading"] = tmp[factor].abs()
        tmp = tmp.sort_values("abs_loading", ascending=False)
        tmp["factor"] = factor

        top_genes_all.append(tmp.head(50))

        safe_factor = sanitize_filename(factor)

        top_pos = tmp.sort_values(factor, ascending=False).head(50)
        top_neg = tmp.sort_values(factor, ascending=True).head(50)
        top_abs = tmp.sort_values("abs_loading", ascending=False).head(100)

        top_pos.to_csv(
            os.path.join(out_dir, f"LDVAE_{exploration_name}_{safe_factor}_top50_positive_loadings.csv"),
            index=False
        )

        top_neg.to_csv(
            os.path.join(out_dir, f"LDVAE_{exploration_name}_{safe_factor}_top50_negative_loadings.csv"),
            index=False
        )

        top_abs.to_csv(
            os.path.join(out_dir, f"LDVAE_{exploration_name}_{safe_factor}_top100_abs_loadings.csv"),
            index=False
        )

    top_genes_ldvae = pd.concat(top_genes_all, axis=0)

    top_genes_ldvae.to_csv(
        os.path.join(out_dir, f"LDVAE_{exploration_name}_top50_genes_per_factor.csv"),
        index=False
    )


def summarize_ldvae_factors(adata, exploration_name, out_dir):
    print("\n\n################################")
    print(f"RELACIÓN FACTORES LDVAE - TEJIDO: {exploration_name}")
    print("################################")

    ldvae_df = pd.DataFrame(
        adata.obsm["X_LDVAE"],
        index=adata.obs_names,
        columns=[f"LDVAE_{i+1}" for i in range(adata.obsm["X_LDVAE"].shape[1])]
    )

    for col in ["region_long", "patient"]:
        ldvae_df[col] = adata.obs[col].values

    ldvae_df.to_csv(
        os.path.join(out_dir, f"LDVAE_{exploration_name}_factors_with_metadata.csv")
    )

    factor_cols = [c for c in ldvae_df.columns if c.startswith("LDVAE_")]

    summary_by_region = ldvae_df.groupby("region_long")[factor_cols].agg(
        ["mean", "median", "std"]
    )

    summary_by_region.to_csv(
        os.path.join(out_dir, f"LDVAE_{exploration_name}_factor_summary_by_region.csv")
    )

    summary_by_patient_region = (
        ldvae_df
        .groupby(["patient", "region_long"])[factor_cols]
        .mean()
        .reset_index()
    )

    summary_by_patient_region.to_csv(
        os.path.join(out_dir, f"LDVAE_{exploration_name}_factor_means_by_patient_region.csv"),
        index=False
    )

    region_order = [
        "Tumor Primario",
        "Metastasis Nodo Linfático"
    ]

    region_order = [
        x for x in region_order
        if x in ldvae_df["region_long"].unique()
    ]

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
            os.path.join(out_dir, f"boxplot_LDVAE_{exploration_name}_{safe_factor}_by_region.png"),
            dpi=300
        )

        plt.close()


# ============================================================
# EJECUCIÓN DE LAS DOS EXPLORACIONES
# ============================================================

for exploration_name, exploration_config in explorations.items():
    print("\n\n============================================================")
    print(f"INICIANDO EXPLORACIÓN: {exploration_name}")
    print(exploration_config["description"])
    print("============================================================")

    out_dir = os.path.join(base_out_dir, exploration_name)
    os.makedirs(out_dir, exist_ok=True)

    try:
        adata_scvi, model_scvi = run_scvi(
            input_h5ad=input_h5ad,
            exploration_name=exploration_name,
            exploration_config=exploration_config,
            out_dir=out_dir
        )
    except Exception as e:
        print(f"\nERROR en scVI para {exploration_name}: {e}")

    try:
        adata_ldvae, model_ldvae = run_ldvae(
            input_h5ad=input_h5ad,
            exploration_name=exploration_name,
            exploration_config=exploration_config,
            out_dir=out_dir
        )
    except Exception as e:
        print(f"\nERROR en LDVAE para {exploration_name}: {e}")


print("\n\n############################################################")
print("DONE: exploraciones no_correction y patient_corrected completadas")
print("############################################################")
