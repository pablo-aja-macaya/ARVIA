from arvia.utils.aeruginosa_snippy import (
    filter_snippy_result,
    filter_synonymous_variants,
    AERUGINOSA_GENES,
)
from arvia.utils.aeruginosa_polymorphisms import PAERUGINOSA_POLYMORPHISMS
from arvia.utils.aeruginosa_truncations import check_truncations
import pandas as pd
import glob
from pathlib import Path
import numpy as np


def concat_files(file_list, selected_bcs=[], header="infer"):
    l = []
    for i in file_list:
        folder = Path(i).parent
        bc = Path(folder).name
        if selected_bcs and bc not in selected_bcs:
            continue
        try:
            temp_df = pd.read_csv(i, sep=None, header=header, engine="python")
            temp_df["bc"] = bc
            l.append(temp_df)

        except Exception as e:
            print(e)
            print("Error: File is not tab or comma delimited, could be empty")

    data = pd.concat(l)
    return data


def apply_paeruginosa_polymorphisms(row, d=PAERUGINOSA_POLYMORPHISMS):
    effect = row["EFFECT"].split(" ")
    prot_effect = [i for i in effect if i.startswith("p.")]
    assert len(prot_effect) <= 1
    if prot_effect:
        prot_effect = prot_effect[0].replace("p.", "")
        if d.get(row["LOCUS_TAG"]) and prot_effect in d.get(row["LOCUS_TAG"]):
            row["polymorphic_substitution"] = "possible polymorphism"

    return row


def get_mutation_props(value):
    try:
        mut_allele_count, ref_allele_count = [int(i.split(":")[1]) for i in value.split(" ")]
        mutation_prop = 100 * mut_allele_count / (mut_allele_count + ref_allele_count)
        return mutation_prop, mut_allele_count
    except Exception as e:
        return None, None


def color_cells(value):
    if value == "-":
        return "text-align: left"
    else:
        mutation_prop, mut_allele_count = get_mutation_props(value)
        if mutation_prop:
            color = None
            if mutation_prop >= 75 and mut_allele_count >= 5:
                color = "#069a2e"
            else:
                color = "#ffaf19"
            return f"background-color : {color}; color: black" + "; text-align: left"
        elif "PMF:" in value:
            color = "#ffaf19"
            return f"background-color : {color}; color: black" + "; text-align: left"
        else:
            return "text-align: left"


def color_cells_v2(value):
    background_color = ""
    s = {"*", "fs", "del", "dup", "ins", "?", "trunc", "possible_missing_feature"}
    search = [True for i in s if i in value]
    if len(search) >= 1:
        background_color = "#ffd17a"
        return f"background-color : {background_color}; color: black" + "; text-align: left"
    elif value.startswith("NZC"):
        d = {i.split("=")[0]: float(i.split("=")[1][:-1]) for i in value.split(" ")}
        if d["NZC"] < 100:
            background_color = "#ffd17a"
            return f"background-color : {background_color}; color: black" + "; text-align: left"
        elif d["NZC"] >= 100 and d["Depth"] < 10:
            background_color = " #ffea69"
            return f"background-color : {background_color}; color: black" + "; text-align: left"
    else:
        return "text-align: left"


def highlight_header(s):
    return [
        "background-color: #2a6099; color: white; font-weight: bold; text-align: center; vertical-align: center; border: 1.2px solid black;"
        for v in s
    ]


def prepare_mutations_for_wide_pivot(data: pd.DataFrame):
    if len(data) >= 1:
        # Normalize positions based on strands (orders mutations)
        data["normalized_pos"] = data["POS"] * (data["STRAND"] + "1").astype(int)
        data = data.sort_values(["CHROM", "LOCUS_TAG", "normalized_pos"], ascending=True).reset_index()

        # Create gene id
        data["gene_id"] = data.apply(
            lambda row: f"{row['LOCUS_TAG']} \n{row['GENE']}"
            if row["GENE"]
            else f"{row['LOCUS_TAG']} \n{row['PRODUCT']}",
            axis=1,
        )

        # Extract mutation parts (nucl, prot...)
        data["mut"] = data["EFFECT"].str.findall("(.*) c\.(.*) p\.(.*)")
        data["mut"] = data.apply(lambda row: row["mut"][0] if row["mut"] else ("", "", row["EFFECT"]), axis=1)
        data[["mutation_type", "mutation_nucl", "mutation_prot"]] = pd.DataFrame(
            data["mut"].to_list(), columns=["mutation_type", "mutation_nucl", "mutation_prot"]
        )

        # Check mutation heterogeneity and depth
        data["mut_qc"] = (
            data["EVIDENCE"]
            .apply(get_mutation_props)
            .apply(
                lambda x: f"(Fails QC: {round(x[0],2)}%, {x[1]}x)"
                if x[0] is not None and (x[0] < 75 or x[1] < 5)
                else ""
            )
        )
        data["mutation_prot"] = (data["mutation_prot"] + " " + data["mut_qc"]).str.strip()

        # Add polymorphisms info
        if "polymorphic_substitution" in data.columns:
            data["mutation_prot"] = data.apply(
                lambda row: f"{row['mutation_prot']} (POLY)"
                if row["polymorphic_substitution"] is not np.nan and row["polymorphic_substitution"] != "-"
                else row["mutation_prot"],
                axis=1,
            )
    else:  # no mutations
        data[
            [
                "normalized_pos",
                "gene_id",
                "mut",
                "mutation_type",
                "mutation_nucl",
                "mutation_prot",
                "mut_qc",
                "mutation_prot",
            ]
        ] = ""

    return data


def get_default_snippy_combination(
    input_files: list, output_file: str = None, selected_bcs: list = [], paeruginosa: bool = False
):

    # -- Concat files --
    df = concat_files(input_files, selected_bcs=selected_bcs)
    if selected_bcs:
        bcs = selected_bcs
        assert len([i for i in selected_bcs if i not in list(df["bc"].unique())]) == 0, "Lacking one input"
    else:
        bcs = list(df["bc"].unique())

    # -- Filter synonymous --
    df = filter_synonymous_variants(df)

    # -- Run paeruginosa filters --
    if paeruginosa:
        df = filter_snippy_result(df).apply(
            apply_paeruginosa_polymorphisms, axis=1
        )  # ATTENTION ONLY FOR PSEUDOMONAS!!!

    # -- Format and copy for later --
    df = df.fillna("-")
    df_copy = df.copy().reset_index()

    # -- Pivot wide --
    index_cols = ["CHROM", "LOCUS_TAG", "POS", "GENE", "PRODUCT", "EFFECT"]
    if "polymorphic_substitution" in df.columns:
        index_cols += ["polymorphic_substitution"]
    pivoted_df = df.pivot(
        index=index_cols,
        columns="bc",
        values="EVIDENCE",
    )

    # -- Paint --
    bcs_with_no_muts = [i for i in bcs if i not in pivoted_df.columns]
    for i in bcs_with_no_muts:
        pivoted_df[i] = "-"

    columns_to_paint = bcs
    pivoted_df = pivoted_df.fillna("-")
    pivoted_df = pivoted_df.style.applymap(color_cells, subset=pd.IndexSlice[:, columns_to_paint])
    pivoted_df = pivoted_df.apply_index(highlight_header, axis="columns", level=[0])

    # Writer object
    if output_file:
        writer = pd.ExcelWriter(output_file, engine="xlsxwriter")

        # Convert the styled dataframe to an XlsxWriter Excel object in specific sheet
        pivoted_df.to_excel(writer, sheet_name="Sheet1")

        # Select sheet and apply formatting
        workbook = writer.book
        worksheet = writer.sheets["Sheet1"]
        # worksheet.autofit()  # autofit row widths
        # worksheet.set_row(1, 45)  # height of row
        # worksheet.set_row(2, 2, [], {"hidden": True})  # row is hidden
        worksheet.set_column(4, 5, 30)  # width of columns 4-5
        worksheet.set_column(7, 7 + len(pivoted_df.columns), 20)  # width of sample columns
        worksheet.freeze_panes(1, 7)  # freeze first row and first 7 column

        # Save
        workbook.close()

    return df_copy, bcs


def paeruginosa_combine_all_mutational_res(
    default_df: pd.DataFrame,
    oprd_fs: list,
    oprd_refs_fs: list,
    gene_coverage_fs: list,
    truncation_fs: list,
    bcs: list,
    output_file: str,
    bcs_without_assembly: list = [],
    filter_poly: bool = False,
):
    if not oprd_fs:
        raise Exception("Expected oprD VC with closest ref. files")

    temp = prepare_mutations_for_wide_pivot(default_df.copy())

    if filter_poly:
        temp = temp[temp["polymorphic_substitution"] != "possible polymorphism"]

    # ---- Prepare oprd closest reference muts ----
    oprd_df = concat_files(oprd_fs)
    oprd_df = filter_synonymous_variants(oprd_df)
    oprd_df = prepare_mutations_for_wide_pivot(oprd_df)
    oprd_df["CHROM"] = oprd_df["CHROM"].str.replace("Pseudomonas_aeruginosa_", "")
    oprd_df = oprd_df[["bc", "CHROM", "gene_id", "mutation_prot"]]

    # -- Pivot all muts --
    temp["gene_id"] = temp["gene_id"] + " \n(Snippy)"
    temp_pivot = (
        temp.drop_duplicates()
        .groupby(["bc", "gene_id"], sort=False, dropna=False)["mutation_prot"]
        .apply(lambda x: ", ".join(list(x)))
        .reset_index(name="muts")
        .pivot(index="bc", columns=["gene_id"], values="muts")
    )
    bcs_with_no_muts = [i for i in bcs if i not in temp_pivot.reset_index()["bc"].to_list()]
    assert len(bcs_with_no_muts) == 0, f"Unexpected samples without mutations: {bcs_with_no_muts}"

    # -- Pivot closest oprD muts --
    temp_oprd = (
        oprd_df.drop_duplicates()
        .groupby(["bc", "gene_id"], sort=False, dropna=False)["mutation_prot"]
        .apply(lambda x: ", ".join(list(x)))
        .reset_index(name="muts")[
            [
                "bc",
                "muts",
            ]
        ]
    )
    # add reference ids on the side because sometimes no mutations are detected in oprd_fs and ref is lost
    oprd_refs = concat_files(oprd_refs_fs, header=None).reset_index(drop=True)
    oprd_refs.columns = ["id", "ref", "CHROM", "identity", "coverage", "bc"]
    temp_oprd = (
        pd.merge(oprd_refs[["bc", "CHROM"]], temp_oprd, on="bc", how="left")
        .rename(
            columns={
                "CHROM": "PA0958-alt \noprD \nclosest reference used",
                "muts": "PA0958-alt \noprD \nclosest reference",
            }
        )
        .fillna("-")
    )

    # -- Truncations based on assembly --
    l = []
    for i in truncation_fs:
        bc = Path(i).parent.name
        temp = check_truncations(i)
        temp["bc"] = bc
        l.append(temp)

    truncations_df = pd.concat(l)
    truncations_df = truncations_df[["bc", "locus_tag", "gene", "comment"]].drop_duplicates()
    truncations_df["locus_gene"] = (
        truncations_df["locus_tag"].astype(str)
        + " \n"
        + truncations_df["gene"].astype(str)
        + " \n"
        + "(Assembly BLAST)"
    )
    temp_truncation = truncations_df.pivot(index="bc", columns="locus_gene", values="comment")

    # ---- Gene coverage ----
    gene_coverage = concat_files(gene_coverage_fs)
    gene_coverage = gene_coverage[["bc", "locus_tag", "non_zero_coverage", "mean_depth"]]
    gene_coverage["comment"] = gene_coverage.apply(
        lambda row: f"NZC={round(100*row['non_zero_coverage'],2)}% Depth={round(row['mean_depth'],2)}x", axis=1
    )
    # gene_coverage["locus_tag"] = gene_coverage["locus_tag"].str.split(".", expand=True)[0]
    gene_coverage = pd.merge(
        pd.DataFrame(AERUGINOSA_GENES.items(), columns=["locus_tag", "gene"]), gene_coverage, on="locus_tag"
    )
    gene_coverage["locus_gene"] = (
        gene_coverage["locus_tag"].astype(str) + " \n" + gene_coverage["gene"].astype(str) + " \n" + "(Gene coverage)"
    )
    temp_gene_cov = gene_coverage.pivot(index="bc", columns="locus_gene", values="comment")

    # -- Merge --
    xxx = pd.merge(temp_pivot, temp_oprd, on="bc", how="left")
    xxx = pd.merge(xxx, temp_truncation, on="bc", how="left")
    xxx = pd.merge(xxx, temp_gene_cov, on="bc", how="left")
    xxx = xxx[["bc"] + [i for i in sorted(xxx.columns) if i != "bc"]]
    xxx[xxx.isna()] = "-"

    # md = pd.read_csv(
    #     "/media/usuario/15d7f455-39ea-40e7-add2-de57c58767eb/analyses/IMIREL/revision_paper_anacandela_imirel/correlation.tsv",
    #     sep="\t",
    # )
    # md = md.rename(columns={"Isolate ID in Database": "bc"})[["CODE", "bc"]]
    # md["bc"] = md["bc"].str.strip()
    # xxx = pd.merge(xxx, md, on="bc", how="left")

    # md = pd.read_excel("/home/usuario/Descargas/ARGA.xlsx")
    # md = md.rename(
    #     columns={
    #         "ID fastq": "bc",
    #         "Profundidad (mediana)": "depth",
    #         "Integridad (%)": "completeness",
    #         "Contaminación (%)": "contamination",
    #     }
    # )[["bc", "ST", "depth", "completeness", "contamination"]]
    # md["bc"] = md["bc"].astype(str)
    # xxx["bc"] = xxx["bc"].astype(str)
    # xxx = pd.merge(xxx, md, on="bc", how="left")

    index_cols = ["bc"]  # , "CODE", "ST", "depth", "completeness", "contamination"
    xxx = xxx[index_cols + [i for i in xxx.columns if i not in index_cols]]

    # xxx.style.applymap(
    #     color_cells_v2,
    #     subset=pd.IndexSlice[
    #         :,
    #         [c for c in xxx.columns if c not in index_cols],
    #     ],
    # ).to_excel("/home/usuario/Proyectos/Results/test7.xlsx")

    yyy = pd.melt(xxx, id_vars=index_cols, var_name="section")
    yyy["locus_tag"] = yyy["section"].str.split(" \n", expand=True)[0]

    # Fill values of samples without blast result with something that indicates the sample did not have an assembly
    if bcs_without_assembly:
        yyy.loc[
            (yyy["section"].str.contains("Assembly BLAST")) 
            & (yyy["bc"].isin(bcs_without_assembly)) 
            & (yyy["value"]=="-"), 
            "value"
        ] = "No assembly given"
    
    # Save this table, which is the long version of the final table and contains everthing, with no " \n"
    zzz = yyy.copy()
    zzz["section"] = zzz["section"].str.replace(" \n", "__").str.replace("(","").str.replace(")","")
    zzz[["bc","locus_tag","section","value"]].to_csv(Path(Path(output_file).parent, "full_long.tsv"), sep="\t", index=None)

    # Pivot and apply style
    temp = (
        yyy.pivot(index=index_cols, columns=["locus_tag", "section"], values="value")
        .style.applymap(
            color_cells_v2,
        )
        .apply_index(highlight_header, axis="columns", level=[0, 1])
        .apply_index(highlight_header, axis="index")
    )

    # Writer object
    writer = pd.ExcelWriter(output_file, engine="xlsxwriter")

    # Convert the styled dataframe to an XlsxWriter Excel object in specific sheet
    temp.to_excel(writer, sheet_name="Sheet1")

    # Select sheet and apply formatting
    workbook = writer.book
    worksheet = writer.sheets["Sheet1"]
    worksheet.autofit()  # autofit row widths
    worksheet.set_row(1, 45)  # height of row
    worksheet.set_row(2, 2, [], {"hidden": True})  # row is hidden # FIXME
    worksheet.set_column(0, 0, 15)  # width of first column
    worksheet.freeze_panes(2, 1)  # freeze first 2 rows and first column

    # Save
    workbook.close()

    return yyy

