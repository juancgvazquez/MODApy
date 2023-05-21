import logging
import multiprocessing as mp
import os
from collections import OrderedDict

from MODApy.cfg import cfg

import cyvcf2

import matplotlib
import matplotlib.pyplot as plt

import matplotlib_venn as venn

import numpy as np

import pandas as pd

matplotlib.use("agg")
logger = logging.getLogger(__name__)


class ParsedVCF(pd.DataFrame):
    """
    A subclass of pandas DataFrame representing parsed VCF data.

    Parameters
    ----------
    pd.DataFrame : pandas DataFrame
        The pandas DataFrame to be subclassed.

    Attributes
    ----------
    _metadata : list of str
        A list containing the names of metadata attributes associated with this object.

    name : str
        The name of the ParsedVCF object.

    Notes
    -----
    The `ParsedVCF` class is a subclass of pandas DataFrame, and thus has all of the
    same functionality and methods as a regular DataFrame. In addition, this class
    provides several methods for parsing VCF files and extracting information from
    them.
    """

    _metadata = ["name"]

    @property
    def _constructor(self):
        """
        Returns the constructor of the current instance.

        Returns
        -------
        type
            The constructor of the current instance, which is always `ParsedVCF`.
        """
        return ParsedVCF

    @classmethod
    def from_vcf(cls, vcf):
        """
        Method that creates a ParsedVCF1 (a DataFrame) from a vcf file
        Parameters
        ----------
        vcf
            Path to the vcf to parse.
        """

        def aminoChange(value: str):
            """
            Given a string `value`, extract the amino acid change from the
            `HGVS.P` format.

            Parameters
            ----------
                value (str): A string representing the `HGVS.P` format.

            Returns
            -------
                str: A string representing the amino acid change or
                "." if not applicable.
            """
            try:
                value = value.replace("p.", "")
                if value[:3] != value[-3:]:
                    return "CHANGE"
                else:
                    return "."
            except Exception:
                return "."

        def divide(x, y):
            """
            Method to divide x on y, needed for dividing freqs.
            Parameters
            ----------
            x
                The dividend
            y
                The divisor
            Returns result or x.
            """
            try:
                return float(x) / y
            except Exception:
                return x

        def parse_vcf_file(vcf):
            """
            Parse a VCF file and return a dictionary of variant information.

            Parameters
            ----------
            vcf : str
                The path to the VCF file to be parsed.

            Returns
            -------
            variants_dict : OrderedDict
                An ordered dictionary of variants with the keys in the format
                CHROM+POS+REF+ALT.
                The values of the dictionary are dictionaries containing information
                about each variant,
                including its ID, QUAL, FILTER, and any additional INFO fields present
                in the VCF file.
            name : str
                The name of the first sample in the VCF file, or the name of the
                VCF file if no samples are present.
            pVCF : cyvcf2.Reader
                A cyvcf2.Reader object representing the VCF file.

            Raises
            ------
            IOError
                If the input VCF file cannot be found or opened.
            """
            logger.info("Parsing VCF File. %s" % vcf)
            pVCF = cyvcf2.Reader(vcf)
            try:
                name = pVCF.samples[0]
            except Exception:
                name = vcf.split("/")[-1]
            variants_dict = OrderedDict()
            for variant in pVCF:
                variants_dict[
                    variant.CHROM
                    + "+"
                    + str(variant.POS)
                    + "+"
                    + variant.REF
                    + "+"
                    + ",".join(variant.ALT)
                ] = {
                    "ID": variant.ID,
                    "QUAL": variant.QUAL,
                    "FILTER": variant.FILTER,
                }
                variants_dict[
                    variant.CHROM
                    + "+"
                    + str(variant.POS)
                    + "+"
                    + variant.REF
                    + "+"
                    + ",".join(variant.ALT)
                ].update({k: v for (k, v) in variant.INFO})
            return variants_dict, name, pVCF

        def create_dataframe(variants_dict):
            """
            Create a pandas DataFrame from the variants dictionary.

            Parameters
            ----------
            variants_dict : dict
                A dictionary of variants data.

            Returns
            -------
            pandas.DataFrame
                A DataFrame containing variants data.

            """
            df1 = pd.DataFrame.from_dict(variants_dict, orient="index")
            del variants_dict
            df1.index = df1.index.str.split("+", expand=True)
            df1.index.names = ["CHROM", "POS", "REF", "ALT"]
            df1.reset_index(inplace=True)
            return df1

        def split_alternate_alleles(df):
            """
            Splits rows with multiple alternate alleles into separate rows.

            Parameters
            ----------
            df : pandas.DataFrame
                The input DataFrame with the genotype data to be processed.

            Returns
            -------
            pandas.DataFrame
                A DataFrame with the same columns as `df`, where rows with multiple
                alternate alleles have been split into separate rows.
            """
            splitdf = df.loc[df["ALT"].str.contains(",")].copy()
            if len(splitdf) > 0:
                ALT = (
                    splitdf["ALT"]
                    .astype(str)
                    .str.split(",", n=1, expand=True)
                    .stack()
                    .rename("ALT")
                )
                ALT.index = ALT.index.droplevel(-1)
                ALT = ALT.to_frame()
                splitdf = splitdf.join(ALT, lsuffix="_x", rsuffix="_y")
                del ALT
                splitdf["ALT"] = splitdf["ALT_y"].combine_first(splitdf["ALT_x"])
                splitdf.drop(columns=["ALT_y", "ALT_x"], inplace=True)
                splitdf.reset_index(inplace=True)
                splitdf.drop(columns="index", inplace=True)
            odd = splitdf.iloc[::2].copy()
            even = splitdf.iloc[1::2].copy()
            splitlist = [
                "ID",
                "AC",
                "AF",
                "SAMPLES_AF",
                "MLEAC",
                "MLEAF",
                "VARTYPE",
                "dbSNPBuildID",
            ]
            splitlist = [x for x in splitlist if x in df.columns]
            splitlist += [x for x in df.columns if x.startswith(("1000", "CLINVAR"))]
            for col in splitlist:
                odd[col] = odd[col].astype(str).str.split(",", n=1).str[0]
                even[col] = even[col].apply(
                    lambda x: x
                    if len(str(x).split(",")) <= 1
                    else str(x).split(",", maxsplit=1)[1]
                )
            splitdf = (
                pd.concat([odd, even])
                .sort_index()
                .replace(to_replace=[r"\(", r"\)"], value="", regex=True)
            )
            del odd, even
            splitdf = splitdf[["CHROM", "POS", "REF", "ALT"] + splitlist]
            df = df.merge(splitdf, on=["CHROM", "POS", "REF"], how="left")
            splitlist.append("ALT")
            xlist = [x + "_x" for x in splitlist]
            ylist = [y + "_y" for y in splitlist]
            del splitdf
            for col in splitlist:
                df[col] = df[col + "_y"].combine_first(df[col + "_x"])
            del splitlist
            df.drop(columns=xlist + ylist, inplace=True)
            del xlist, ylist
            df["POS"] = df["POS"].astype(int)
            return df

        def handle_annotations(df, pVCF):
            """
            Parses the 'ANN' column in a pandas DataFrame and extracts functional
            annotations as separate columns.

            Parameters:
            -----------
            df : pandas.DataFrame
                Input DataFrame with 'ANN' column containing functional annotations
                separated by commas.
            pVCF : cyvcf2.Reader
                A cyvcf2.Reader object representing the VCF file.

            Returns:
            --------
            pandas.DataFrame
                A new DataFrame with functional annotations as separate columns,
                joined with the original DataFrame.
            """
            if "ANN" in df.columns:
                anndf = df["ANN"]
                annhead = pVCF.get_header_type("ANN")["Description"].strip(
                    '"Functional annotations: \'"'
                )
                annheaderlist = [x.strip() for x in annhead.split("|")]
                anndf = anndf.str.split(",", expand=True).stack()
                anndf = anndf.str.split("|", expand=True)
                anndf.columns = annheaderlist
                df.drop(columns="ANN", inplace=True)
                anndf.index = anndf.index.droplevel(1)
                df = df.join(anndf, how="inner")
                del anndf
                del annhead
                del annheaderlist
            return df

        def prioritize_variants(df, IMPACT_SEVERITY=None):
            """
            Sort variants in a pandas DataFrame according to their severity.

            Parameters
            ----------
            df : pandas.DataFrame
                The DataFrame containing the variants to be sorted.
            IMPACT_SEVERITY : dict, optional
                A dictionary that maps each variant type to its severity score.
                The keys of the dictionary should be the names of the variant types
                (e.g., 'missense_variant'), and the values should be integers
                representing the severity score. If not provided, the default values for
                severity scores will be used.

            Returns
            -------
            pandas.DataFrame
                A new DataFrame with the same columns as the input DataFrame, but with
                the variants sorted by their severity.

            Notes
            -----
            This function assumes that the input DataFrame has columns named 'CHROM',
            'POS', 'REF', 'ALT', 'Annotation', and 'HGVS.c'. The 'Annotation' column
            should contain the variant types
            (e.g., 'missense_variant&splice_region_variant'), and the 'HGVS.c' column
            should contain the HGVS coding sequence notation for each variant
            (e.g., 'NM_001005353.2:c.43A>G'). Variants with a null HGVS notation will be
            sorted to the end.
            """
            if IMPACT_SEVERITY is None:
                IMPACT_SEVERITY = {
                    "exon_loss_variant": 1,
                    "frameshift_variant": 2,
                    "stop_gained": 3,
                    "stop_lost": 4,
                    "start_lost": 5,
                    "splice_acceptor_variant": 6,
                    "splice_donor_variant": 7,
                    "disruptive_inframe_deletion": 8,
                    "inframe_insertion": 9,
                    "disruptive_inframe_insertion": 10,
                    "inframe_deletion": 11,
                    "missense_variant": 12,
                    "splice_region_variant": 13,
                    "stop_retained_variant": 14,
                    "initiator_codon_variant": 15,
                    "synonymous_variant": 16,
                    "start_retained": 17,
                    "coding_sequence_variant": 18,
                    "5_prime_UTR_variant": 19,
                    "3_prime_UTR_variant": 20,
                    "5_prime_UTR_premature_start_codon_gain_variant": 21,
                    "intron_variant": 22,
                    "non_coding_exon_variant": 23,
                    "upstream_gene_variant": 24,
                    "downstream_gene_variant": 25,
                    "TF_binding_site_variant": 26,
                    "regulatory_region_variant": 27,
                    "intergenic_region": 28,
                    "transcript": 29,
                }
            if 'Annotation' in df.columns:
                df["sorter"] = (
                    df["Annotation"].str.split("&").str[0].replace(IMPACT_SEVERITY)
                )
                df.loc[df["HGVS.c"].str.contains("null"), "HGVS.c"] = None
                df["sorter2"] = [x[0] == x[1] for x in zip(df["ALT"], df["Allele"])]
                df = df.sort_values(
                    by=["CHROM", "POS", "sorter2", "sorter"],
                    ascending=[True, True, False, True],
                ).drop_duplicates(["CHROM", "POS", "REF", "ALT"])
                df.drop(columns=["sorter", "sorter2"], inplace=True)
            return df

        def format_ann_columns(df, pVCF):
            """
            Formats the columns of a pandas DataFrame containing variant annotation
            data.

            Parameters
            ----------
            df : ParsedVCF (pandas.DataFrame extension)
                A DataFrame containing variant annotation data.
            pVCF : cyvcf2.Reader
                A cyvcf2.Reader object representing the VCF file.

            Returns
            -------
            ParsedVCF (pandas.DataFrame extension)
                The input DataFrame with formatted columns.

            Notes
            -----
            The function applies the following transformations to the input DataFrame:

            - All column names are converted to uppercase.
            - If the DataFrame contains a column named 'HGVS.P', a new column named
            'AMINOCHANGE'
            is added, which contains the result of calling the `aminoChange` function
            on the 'HGVS.P'
            column.
            - If the DataFrame contains a column named 'HOM', its values are converted
            from boolean to
            categorical ('HOM' and 'HET'). The 'HET' column is dropped if present. The
            column name is changed to 'ZIGOSITY'.
            - If the DataFrame contains a column named 'ESP6500_MAF', new columns named
            'ESP6500_MAF_EA', 'ESP6500_MAF_AA', and 'ESP6500_MAF_ALL' are added,
            containing the values of the corresponding fields in the 'ESP6500_MAF'
            column. The values are converted from strings to floats and divided by 100.
            The 'ESP6500_MAF' column is dropped.
            - If the DataFrame contains a column named 'ESP6500_PH', new columns named
            'POLYPHEN_PRED' and 'POLYPHEN_SCORE' are added, containing the values of the
            corresponding fields in the 'ESP6500_PH' column. The 'POLYPHEN_PRED' values
            are cleaned up by removing trailing dots and commas. The 'POLYPHEN_SCORE'
            values are split on commas and the first element is kept. The 'ESP6500_PH'
            column is dropped.
            - The columns named 'ANNOTATION', 'ANNOTATION_IMPACT', and 'ID' are renamed
            to 'EFFECT', 'IMPACT', and 'RSID', respectively.
            - Columns with numeric data (according to the VCF header) are converted to
            floats or integers, as appropriate. The columns named 'ESP6500_MAF_EA',
            'ESP6500_MAF_AA', and 'ESP6500_MAF_ALL' are also converted to floats.
            - The DataFrame is rounded to 6 decimal places.
            - If the DataFrame contains a column named 'CLINVAR_CLNSIG', its values are
            replaced with their corresponding meanings according to the
            `clinvartranslation` dictionary.

            """
            df.columns = df.columns.str.upper()
            if "HGVS.P" in df.columns:
                df["AMINOCHANGE"] = df["HGVS.P"].apply(aminoChange)
            if "HOM" in df.columns:
                df["HOM"] = df["HOM"].replace({True: "HOM", np.nan: "HET", None: "HET"})
                df.drop(columns="HET", inplace=True, errors="ignore")
                df.rename(columns={"HOM": "ZIGOSITY"}, inplace=True)
            if "ESP6500_MAF" in df.columns:
                df[["ESP6500_MAF_EA", "ESP6500_MAF_AA", "ESP6500_MAF_ALL"]] = df[
                    "ESP6500_MAF"
                ].str.split(",", expand=True)
                df["ESP6500_MAF_EA"] = df["ESP6500_MAF_EA"].apply(divide, args=(100,))
                df["ESP6500_MAF_AA"] = df["ESP6500_MAF_AA"].apply(divide, args=(100,))
                df["ESP6500_MAF_ALL"] = df["ESP6500_MAF_ALL"].apply(divide, args=(100,))
                df.drop(columns=["ESP6500_MAF"], inplace=True)
            if "ESP6500_PH" in df.columns:
                df[["POLYPHEN_PRED", "POLYPHEN_SCORE"]] = df["ESP6500_PH"].str.split(
                    ":", n=1, expand=True
                )
                df["POLYPHEN_PRED"] = df["POLYPHEN_PRED"].str.strip(".").str.strip(".,")
                df["POLYPHEN_SCORE"] = df["POLYPHEN_SCORE"].str.split(",").str[0]
                df.drop(columns=["ESP6500_PH"], inplace=True)
            df.rename(
                columns={
                    "ANNOTATION": "EFFECT",
                    "ANNOTATION_IMPACT": "IMPACT",
                    "ID": "RSID",
                },
                inplace=True,
                errors="ignore",
            )
            numcols = list()
            for x in pVCF.header_iter():
                if x.type == "INFO":
                    if x["Type"] in ["Float", "Integer"]:
                        numcols.append(x["ID"])
            numcols += ["ESP6500_MAF_EA", "ESP6500_MAF_AA", "ESP6500_MAF_ALL"]
            numcols = list(
                set([x.upper() for x in numcols for y in df.columns if x.upper() == y])
            )
            df[numcols] = df[numcols].apply(pd.to_numeric, errors="coerce", axis=1)
            df = df.round(6)

            if "CLINVAR_CLNSIG" in df.columns:
                clinvartranslation = {
                    "255": "other",
                    "0": "Uncertain significance",
                    "1": "not provided",
                    "2": "Benign",
                    "3": "Likely Benign",
                    "4": "Likely pathogenic",
                    "5": "Pathogenic",
                    "6": "drug response",
                    "7": "histocompatibility",
                }
                for k, v in clinvartranslation.items():
                    df["CLINVAR_CLNSIG"] = df["CLINVAR_CLNSIG"].str.replace(k, v)
            return df

        def clean_df(df):
            """
            Replace missing and empty values in the DataFrame with the '.' character.
            Convert the 'POS' column to integer type.

            Parameters
            ----------
            df : pandas.DataFrame
                Input DataFrame to be cleaned.

            Returns
            -------
            pandas.DataFrame
                Cleaned DataFrame with replaced missing and empty values and 'POS'
                column as integer.

            """
            df.replace(["nan", "", np.nan], ".", inplace=True)
            df.replace(to_replace=[None], value=".", inplace=True, regex=True)
            df = df.astype("str")
            df["POS"] = df["POS"].astype(int)
            return df

        variants_dict, name, pVCF = parse_vcf_file(vcf)
        df1 = create_dataframe(variants_dict)
        del variants_dict
        df1 = split_alternate_alleles(df1)
        df1 = handle_annotations(df1, pVCF)
        df1 = prioritize_variants(df1)
        df1 = format_ann_columns(df1, pVCF)
        df1 = clean_df(df1)
        df1 = df1.pipe(ParsedVCF)
        df1.name = name
        return df1

    @classmethod
    def mp_parser(cls, *vcfs, cores=int(cfg["GENERAL"]["cores"])):
        """
        Parses multiple VCF files concurrently using multiprocessing.

        Parameters:
        -----------
        `*vcfs` : str
            Paths to the VCF files to be parsed.
        cores : int, optional
            Number of cores to be used for multiprocessing. If not specified,
            the number of cores specified in the configuration file will be used.

        Returns:
        --------
        list of ParsedVCF
            A list of ParsedVCF objects parsed from the input VCF files.
        """
        if len(vcfs) < 1:
            logger.error("No vcfs provided!")
            exit(1)
        elif len(vcfs) == 1:
            pvcfs = list()
            pvcfs.append(ParsedVCF.from_vcf(vcfs[0]))
        else:
            try:
                [x + "" for x in vcfs]
            except Exception:
                logger.error("All mp_parser args must be strings")
            else:
                logger.info("Starting Multi-Parser")
                if cores is None:
                    cores = mp.cpu_count()
                if cores > 1:
                    if len(vcfs) <= cores - 1:
                        pool = mp.Pool(processes=len(vcfs))
                    else:
                        pool = mp.Pool(processes=cores - 1)
                else:
                    pool = mp.Pool(processes=cores)
                pvcfs = pool.map(cls.from_vcf, (x for x in vcfs))
                pool.close()
                pool.join()
        return pvcfs

    def to_macrogen_xls(self, outpath):
        """
        Writes a sorted and filtered DataFrame to an Excel file with columns in the
        format required by Macrogen.

        Parameters:
        -----------
        outpath : str
            The output file path for the Excel file.

        Returns:
        --------
        None
        """
        macrogen_cols = [
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "DP",
            "AD",
            "QUAL",
            "MQ",
            "Zygosity",
            "FILTER",
            "Effect",
            "IMPACT",
            "Gene_Name",
            "Feature_Type",
            "Feature_ID",
            "Transcript_BioType",
            "Rank/Total",
            "HGVS.c",
            "HGVS.p",
            "REF_AA",
            "ALT_AA",
            "cDNA_pos/cDNA_length",
            "CDS_pos/CDS_length",
            "AA_pos/AA_length",
            "Distance",
            "dbSNP142_ID",
            "1000Gp3_AF",
            "1000Gp3_AFR_AF",
            "1000Gp3_AMR_AF",
            "1000Gp3_EAS_AF",
            "1000Gp3_EUR_AF",
            "1000Gp3_SAS_AF",
            "ESP6500_MAF_EA",
            "ESP6500_MAF_AA",
            "ESP6500_MAF_ALL",
            "SIFT_score",
            "SIFT_pred",
            "Polyphen2_HDIV_score",
            "Polyphen2_HDIV_pred",
            "Polyphen2_HVAR_score",
            "Polyphen2_HVAR_pred",
            "CLINVAR_CLNSIG",
            "CLINVAR_CLNDSDB",
            "CLINVAR_CLNDSDBID",
            "CLINVAR_CLNDBN",
            "CLINVAR_CLNREVSTAT",
            "CLINVAR_CLNACC",
        ]
        self["chrsort"] = self["CHROM"].replace({"X": 30, "Y": 40})
        df1 = self.sort_values(["chrsort", "POS"])[
            [x.upper() for x in macrogen_cols if x.upper() in self.columns]
        ].copy()
        df1.to_excel(outpath)

    def panel(self, panel):
        """
        Extracts a subset of variants from the current ParsedVCF object based on a
        panel of genes.

        Parameters
        ----------
        panel : str or path-like object
            Path to an Excel file containing a sheet named "GeneList" with a column
            "GeneSymbol" containing the list of genes to extract.

        Returns
        -------
        panel_df : ParsedVCF object
            A new ParsedVCF object containing only variants from the genes in the panel.

        Raises
        ------
        Exception
            If the "GeneList" sheet cannot be parsed from the panel Excel file.

        Notes
        -----
        - The resulting `panel_df` object will have the same name as the original \
        object.
        - If no variants are found for any gene in the panel, an error message will \
        be logged.
        """
        logger.info("Analyzing Panel")
        try:
            pldf = pd.read_excel(panel, sheet_name="GeneList")
            geneSymbolList = list(pldf.GeneSymbol.unique())
        except Exception:
            logger.error("There was an error parsing GeneList")
            logger.debug("", exc_info=True)
            exit(1)
        panel_df = pd.DataFrame()
        for gene in geneSymbolList:
            panel_df = panel_df.append(self.loc[self["GENE_NAME"] == gene])
        panel_df = panel_df.pipe(ParsedVCF)
        panel_df.name = self.name
        if len(panel_df) < 1:
            logger.error("After running Panel, resulting Dataframe holds no variants")
        return panel_df

    def duos(self, vcf2, VENNPLACE=None):
        """
        Compare two VCF files using CHROM, POS, REF, and ALT columns as the index.

        Parameters
        ----------
        vcf2 : str or ParsedVCF
            The VCF file or ParsedVCF object to compare to.
        VENNPLACE : str, optional
            The path to save the Venn diagram. Defaults to None.

        Returns
        -------
        pandas.DataFrame
            A new DataFrame containing a 'DUOS' column that indicates which file the
            variant belongs to.

        Raises
        ------
        SystemExit
            If both `self` and `vcf2` have a 'VENN' column (meaning they are both duos),
            then the function will raise a `SystemExit`.

        Notes
        -----
        This method checks if the second argument (`vcf2`) is a path to a VCF file or a
        `ParsedVCF` object. It then drops specific columns (`indcols`) from `self` and
        `vcf2`.

        If both `self` and `vcf2` have a 'VENN' column, this means that both files are
        duos. The function cannot currently combine more than one duo, so it will raise
        a `SystemExit`. If `self` has a 'VENN' column, then `vcf2` will be the left
        dataframe. If `vcf2` has a 'VENN' column, then `self` will be the left
        dataframe. Otherwise, `self` will be the left dataframe and `vcf2` will be the
        right dataframe.

        The function then performs an outer merge on the `CHROM`, `POS`, `REF`, and
        `ALT` columns, using the `indicator` parameter to indicate whether a variant
        is in `self`, `vcf2`, or both. It then drops specific columns (`difcols` and
        `eqcols`) and returns the merged dataframe.

        If the `VENNPLACE` parameter is provided, the function will save a Venn
        diagram to that path. The Venn diagram will indicate whether a variant is in
        `self`, `vcf2`, or both.
        """

        # chequeo si el segundo es un vcf parseado o una ruta a un vcf
        def _trios_stats(self, names):
            logger.info("Calculating Trios statistics")
            trios = self.groupby("VENN", sort=False).size()
            A, B, C = names.split(":")
            names = [A, B]
            names.append(A + ":" + B)
            names.append(C)
            names.append(A + ":" + C)
            names.append(B + ":" + C)
            names.append(":".join([A, B, C]))
            trios = trios.reindex(names).fillna(0)
            trios = trios.astype(int)
            triosgraph = plt.figure()
            venn.venn3(trios, set_labels=[A, B, C], set_colors=["b", "r", "g"])
            triosgraph.savefig("./venn.png", dpi=triosgraph.dpi)
            triosgraph.clf()

        def _duos_stats(self, names):
            logger.info("Calculating Duos statistics")
            duos = self.groupby("VENN", sort=False).size()
            A, B = names.split(":")
            names = [A, B, names]
            duos = duos.reindex(names).fillna(0)
            duosgraph = plt.figure()
            venn.venn2(duos, set_labels=[A, B], set_colors=["b", "r"])
            duosgraph.savefig("./venn.png", dpi=duosgraph.dpi)
            duosgraph.clf()

        if isinstance(vcf2, str):
            pvcf2 = ParsedVCF.from_vcf(vcf2)
        elif isinstance(vcf2, ParsedVCF):
            pvcf2 = vcf2

        indcols = [
            "QUAL",
            "FILTER",
            "DP",
            "FS",
            "MQ",
            "SOR",
            "QD",
            "SET",
            "BASEQRANKSUM",
            "CLIPPINGRANKSUM",
            "MQRANKSUM",
            "READPOSRANKSUM",
            "AC",
            "SAMPLES_AF",
            "MLEAC",
            "MLEAF",
            "DBSNPBUILDID",
        ]
        indself = [x for x in indcols if x in self.columns]
        indpvcf2 = [x for x in indcols if x in pvcf2.columns]
        self.drop(columns=indself, inplace=True)
        pvcf2.drop(columns=indpvcf2, inplace=True)
        del indcols, indself, indpvcf2
        # chequeo si alguno es un duos y dropeo columnas individuales
        if ("VENN" in self.columns) & ("VENN" in pvcf2.columns):
            logger.error(
                "Only one argument can be a duos. Function cannot \
                currently combine more than one duos (resulting"
                "in a trios). Both arguments provided where duos or trios"
            )
            exit(1)
        elif "VENN" in self.columns:
            logger.info(
                "Running TRIOS analysis on %s" % ":".join([self.name, pvcf2.name])
            )
            indicator = "TRIOS"
            left = pvcf2
            right = self
        elif "VENN" in pvcf2.columns:
            logger.info(
                "Running TRIOS analysis on %s" % ":".join([self.name, pvcf2.name])
            )
            indicator = "TRIOS"
            left = self
            right = pvcf2
        else:
            logger.info(
                "Running DUOS analysis on %s" % ":".join([self.name, pvcf2.name])
            )
            indicator = "DUOS"
            left = self
            right = pvcf2

        # Hago el merge
        mergedVCF = left.merge(
            right,
            on=["CHROM", "POS", "REF", "ALT"],
            how="outer",
            suffixes=("_" + self.name, "_" + pvcf2.name),
            indicator=indicator,
        )

        # columnas que deberían ser iguales y columnas que podrían ser distintas
        difcols = [
            x.replace("_" + self.name, "")
            for x in mergedVCF.columns
            if "_" + self.name in x
        ]
        eqcols = [
            x
            for x in difcols
            if any(
                y in x for y in ["1000GP3", "CLINVAR", "ESP6500", "RSID", "POLYPHEN"]
            )
        ]
        difcols = set(difcols) - set(eqcols)
        # columnas a eliminar
        dropcols = [x + "_" + self.name for x in difcols]
        dropcols += [x + "_" + pvcf2.name for x in difcols]
        dropcols += [x + "_" + pvcf2.name for x in eqcols]
        dropcols += [x + "_" + self.name for x in eqcols]
        # armo un dataframe para combinar las columnas
        tmp = mergedVCF.loc[mergedVCF[indicator] == "both"][dropcols]
        # combino las que deberían ser iguales
        for col in eqcols:
            mergedVCF[col] = mergedVCF[col + "_" + self.name].combine_first(
                mergedVCF[col + "_" + pvcf2.name]
            )
        # combino las que podrían ser diferentes, si son diferentes (para
        # variantes en ambos archivos, no las combino.
        for col in difcols:
            if all(tmp[col + "_" + self.name] == tmp[col + "_" + pvcf2.name]):
                mergedVCF[col] = mergedVCF[col + "_" + self.name].combine_first(
                    mergedVCF[col + "_" + pvcf2.name]
                )
            else:
                dropcols.remove(col + "_" + self.name)
                dropcols.remove(col + "_" + pvcf2.name)
        # dropeo las columans que combine
        mergedVCF.drop(columns=dropcols, inplace=True)
        # armo la columna indicadora para Duos y Trios
        if indicator == "DUOS":
            mergedVCF["DUOS"].replace(
                {
                    "left_only": self.name,
                    "right_only": pvcf2.name,
                    "both": self.name + ":" + pvcf2.name,
                },
                inplace=True,
            )
            mergedVCF.rename(columns={"DUOS": "VENN"}, inplace=True)
            names = ":".join([self.name, pvcf2.name])
            _duos_stats(mergedVCF, names)
            if VENNPLACE is not None:
                print(VENNPLACE)
                if VENNPLACE == "A":
                    mergedVCF = mergedVCF[mergedVCF["VENN"] == self.name]
                elif VENNPLACE == "B":
                    mergedVCF = mergedVCF[mergedVCF["VENN"] == pvcf2.name]
                elif VENNPLACE == "A:B":
                    mergedVCF = mergedVCF[
                        mergedVCF["VENN"] == self.name + ":" + pvcf2.name
                    ]
                else:
                    logger.error(
                        "VENNPLACE can only be A, B or A:B for a trios analysis"
                    )
                    logger.debug("", exc_info=True)
                    exit(1)
        if indicator == "TRIOS":
            mergedVCF["PATIENT"] = None
            mergedVCF["PATIENT"] = np.where(
                mergedVCF["TRIOS"] == "left_only",
                left.name,
                mergedVCF["PATIENT"],
            )
            mergedVCF["PATIENT"] = np.where(
                (mergedVCF["TRIOS"] == "both"),
                mergedVCF["VENN"] + ":" + left.name,
                mergedVCF["PATIENT"],
            )
            mergedVCF["PATIENT"] = np.where(
                (mergedVCF["TRIOS"] == "right_only"),
                mergedVCF["VENN"],
                mergedVCF["PATIENT"],
            )
            mergedVCF.drop(columns=["VENN", "TRIOS"], inplace=True)
            mergedVCF.rename(
                columns={
                    "PATIENT": "VENN",
                    "ZIGOSITY": "ZIGOSITY_" + left.name,
                },
                inplace=True,
            )
            names = self.name + ":" + pvcf2.name
            _trios_stats(mergedVCF, names)

            if VENNPLACE is not None:
                names = names.split(":")
                if VENNPLACE == "A":
                    mergedVCF = mergedVCF[
                        (mergedVCF["VENN"].str.contains(names[0]))
                        & ~(mergedVCF["VENN"].str.contains(names[1]))
                        & ~(mergedVCF["VENN"].str.contains(names[2]))
                    ]
                elif VENNPLACE == "B":
                    mergedVCF = mergedVCF[
                        ~(mergedVCF["VENN"].str.contains(names[0]))
                        & (mergedVCF["VENN"].str.contains(names[1]))
                        & ~(mergedVCF["VENN"].str.contains(names[2]))
                    ]
                elif VENNPLACE == "C":
                    mergedVCF = mergedVCF[
                        ~(mergedVCF["VENN"].str.contains(names[0]))
                        & ~(mergedVCF["VENN"].str.contains(names[1]))
                        & (mergedVCF["VENN"].str.contains(names[2]))
                    ]
                elif VENNPLACE == "A:B":
                    mergedVCF = mergedVCF[
                        (mergedVCF["VENN"].str.contains(names[0]))
                        & (mergedVCF["VENN"].str.contains(names[1]))
                        & ~(mergedVCF["VENN"].str.contains(names[2]))
                    ]
                elif VENNPLACE == "A:C":
                    mergedVCF = mergedVCF[
                        (mergedVCF["VENN"].str.contains(names[0]))
                        & ~(mergedVCF["VENN"].str.contains(names[1]))
                        & (mergedVCF["VENN"].str.contains(names[2]))
                    ]
                elif VENNPLACE == "B:C":
                    mergedVCF = mergedVCF[
                        ~(mergedVCF["VENN"].str.contains(names[0]))
                        & (mergedVCF["VENN"].str.contains(names[1]))
                        & (mergedVCF["VENN"].str.contains(names[2]))
                    ]
                elif VENNPLACE == "A:B:C":
                    mergedVCF = mergedVCF[
                        (mergedVCF["VENN"].str.contains(names[0]))
                        & (mergedVCF["VENN"].str.contains(names[1]))
                        & (mergedVCF["VENN"].str.contains(names[2]))
                    ]
                else:
                    logger.error(
                        "VENNPLACE can only be A, B, C, A:B, A:C, B:C or A:B:C \
                            in a trios analysis"
                    )
                    logger.debug("", exc_info=True)
                    exit(1)
        if len(mergedVCF) < 1:
            logger.error(
                "After running Duos/Trios, resulting Dataframe does not \
                    hold any variants"
            )
            exit(1)
        mergedVCF.fillna(".", inplace=True)
        mergedVCF = mergedVCF.pipe(ParsedVCF)
        mergedVCF.name = ":".join([self.name, pvcf2.name])
        return mergedVCF

    def general_stats(self):
        """
        Calculates general statistics of the VCF file.

        Returns:
        --------
        vcfstats : pandas DataFrame
            A DataFrame containing the calculated statistics grouped by CHROM,
            VARTYPE, IMPACT, and EFFECT (if "VENN" is in self.columns), or by
            CHROM, ZIGOSITY, VARTYPE, IMPACT, and EFFECT (otherwise). The count
            column shows the number of variants in each group.

        Raises:
        -------
        None

        Example:
        --------
        >>> vcf = VCF("example.vcf")
        >>> stats = vcf.general_stats()
        """
        if "VENN" in self.columns:
            colstats = ["CHROM", "VARTYPE", "IMPACT", "EFFECT"]
        else:
            colstats = ["CHROM", "ZIGOSITY", "VARTYPE", "IMPACT", "EFFECT"]
        if set(colstats).issubset(self.columns):
            logger.info("Calculating General Statistics")
            vcfstats = self.groupby(colstats).size().to_frame(name="count")
            vcfstats = vcfstats.rename_axis(colstats).reset_index()
            vcfstats.name = "stats"
            plt.pie(
                list(vcfstats.groupby("CHROM").size().values),
                labels=vcfstats.groupby("CHROM").size().index.values,
            )
            my_circle = plt.Circle((0, 0), 0.7, color="white")
            chromVars = plt.gcf()
            chromVars.gca().add_artist(my_circle)
            chromVars.savefig("./general.png")
            return vcfstats

    def vcf_to_excel(self, outpath):
        """
        Convert the variant call format (VCF) data in the pandas DataFrame
        to an Excel file format using pandas ExcelWriter.

        Parameters
        ----------
        self : pandas.DataFrame
            A pandas DataFrame containing the VCF data.
        outpath : str
            The path where the resulting Excel file will be written.

        Returns
        -------
        None

        Raises
        ------
        OSError
            If the directory specified by outpath.rsplit('/', maxsplit=1)[0] does not
            exist and cannot be created.

        Notes
        -----
        This method performs the following operations:
        - Creates the directory specified by outpath.rsplit('/', maxsplit=1)[0] if it
        does not exist.
        - Creates a pandas ExcelWriter object with the file specified by outpath.
        - Adds an empty 'VARSOME' column to the DataFrame.
        - Selects the columns in the DataFrame specified by the configuration file's
        OUTPUT.columnsorder setting.
        - Adds columns with 'ZIGOSITY' in their name to the selected columns if 'VENN'
        is in the DataFrame columns and
        'ZIGOSITY' is in the selected columns.
        - Selects only the columns that exist in the DataFrame and the selected columns
        list.
        - Sorts the DataFrame by the first selected column.
        - Creates a new worksheet named 'DATA' in the Excel file.
        - Sets the number format to '0.00000' for all columns in the 'DATA' worksheet.
        - Sets the number format to '###,###,###' for the 'POS' column.
        - Sets the number format to 'General' for the 'RSID' column.
        - Applies conditional formatting to the 'IMPACT' column using four different
        formats based on the text content.
        - Writes the sorted and formatted DataFrame to the 'DATA' worksheet.
        - Redirects the 'RSID' and 'GENE_NAME' columns to URLs if the index of the reset
        DataFrame is less than 32150.
        - Inserts an autofilter to the 'DATA' worksheet.
        - Calculates general statistics for the DataFrame.
        - Creates a new worksheet named 'STATISTICS' in the Excel file.
        - Writes the statistics DataFrame to the 'STATISTICS' worksheet.
        - Inserts an image of the general statistics to the 'STATISTICS' worksheet.
        - Inserts an image of a Venn diagram to the 'STATISTICS' worksheet if the file
        ./venn.png exists.
        """
        os.makedirs(outpath.rsplit("/", maxsplit=1)[0], exist_ok=True)
        output = pd.ExcelWriter(outpath)
        self["VARSOME"] = ""
        cols_selected = cfg["OUTPUT"]["columnsorder"].replace(",", " ").split()
        if "VENN" in self.columns:
            if "ZIGOSITY" in cols_selected:
                cols_selected += [x for x in self.columns if "ZIGOSITY" in x]
        finalcols = [x for x in cols_selected if x in self.columns]
        self = self[finalcols]
        self = self.sort_values(by=finalcols[0])
        workbook = output.book
        datasheet = workbook.add_worksheet("DATA")
        statsheet = workbook.add_worksheet("STATISTICS")

        output.sheets["DATA"] = datasheet

        formatnum = workbook.add_format({"num_format": "0.00000"})
        # for i, col in enumerate(self.columns):
        datasheet.set_column(0, len(self.columns), 15, formatnum)

        formatpos = workbook.add_format({"num_format": "###,###,###"})
        datasheet.set_column(
            finalcols.index("POS"), finalcols.index("POS"), 15, formatpos
        )
        datasheet.set_column(finalcols.index("RSID"), finalcols.index("RSID"), 15)
        # Light red fill with dark red text.
        highformat = workbook.add_format(
            {"bg_color": "#FFC7CE", "font_color": "#9C0006", "bold": True}
        )
        # Light yellow fill with dark yellow text.
        modformat = workbook.add_format(
            {"bg_color": "#FFFF99", "font_color": "#9C6500", "bold": True}
        )
        # Light orange fill with dark orange text.
        moderformat = workbook.add_format(
            {"bg_color": "#FFCC99", "font_color": "#FF6600", "bold": True}
        )
        # Green fill with dark green text.
        lowformat = workbook.add_format(
            {"bg_color": "#C6EFCE", "font_color": "#006100", "bold": True}
        )
        datasheet.conditional_format(
            0,
            cols_selected.index("IMPACT"),
            len(self),
            cols_selected.index("IMPACT"),
            {
                "type": "text",
                "criteria": "containing",
                "value": "HIGH",
                "format": highformat,
            },
        )
        datasheet.conditional_format(
            0,
            cols_selected.index("IMPACT"),
            len(self),
            cols_selected.index("IMPACT"),
            {
                "type": "text",
                "criteria": "containing",
                "value": "MODIFIER",
                "format": modformat,
            },
        )
        datasheet.conditional_format(
            0,
            cols_selected.index("IMPACT"),
            len(self),
            cols_selected.index("IMPACT"),
            {
                "type": "text",
                "criteria": "containing",
                "value": "MODERATE",
                "format": moderformat,
            },
        )
        datasheet.conditional_format(
            0,
            cols_selected.index("IMPACT"),
            len(self),
            cols_selected.index("IMPACT"),
            {
                "type": "text",
                "criteria": "containing",
                "value": "LOW",
                "format": lowformat,
            },
        )
        logger.info("Writing Excel File")
        self.to_excel(
            output,
            sheet_name="DATA",
            merge_cells=False,
            index=False,
            header=True,
        )

        if self.reset_index().index.max() < 32150:
            logger.info("Redirecting IDs and GENEs to URLs")
            try:
                colid = cols_selected.index("RSID")
                colgen = cols_selected.index("GENE_NAME")
                row = 2
                for x in zip(self["RSID"], self["GENE_NAME"]):
                    if type(x[0]) == str:
                        urlrs = "https://varsome.com/variant/hg19/%s"
                        rsvalue = (x[0].replace(";", ",").split(","))[0]
                        datasheet.write_url(
                            "%s%i" % (chr(colid + 65), (row)),
                            urlrs % rsvalue,
                            string=rsvalue,
                        )
                    if type(x[1]) == str:
                        urlgen = "https://www.ncbi.nlm.nih.gov/omim/?term=%s"
                        datasheet.write_url(
                            "%s%i" % (chr(colgen + 65), (row)),
                            urlgen % x[1],
                            string=x[1],
                        )
                    row += 1
            except Exception as e:
                logger.error(e, exc_info=True)
        datasheet.autofilter(0, 0, len(self), len(finalcols))
        stats = self.general_stats()
        output.sheets["STATISTICS"] = statsheet
        try:
            stats.to_excel(output, sheet_name="STATISTICS")
        except Exception as e:
            logger.error(
                "Could not print statistics. Error was {}".format(e),
                exc_info=True,
            )
        try:
            statsheet.insert_image("H2", "./general.png")
        except Exception as e:
            logger.error(
                "Could not print stats graphs. Error was {}".format(e),
                exc_info=True,
            )
        if os.path.isfile("./venn.png"):
            statsheet.insert_image("H25", "./venn.png")
        output.save()
        try:
            os.remove("./general.png")
        except Exception:
            logger.debug("Could not remove general.png")
        try:
            os.remove("./venn.png")
        except Exception:
            logger.debug("Could not remove venn.png")
