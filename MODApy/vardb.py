from MODApy.cfg import configuration

import duckdb
import logging

logger = logging.getLogger(__name__)


class VarDB:
    @property
    def _constructor(self):
        return VarDB

    def __init__(self, conn_type='file', conn_str=configuration.variantsDBPath):
        if conn_type.lower() == "file":
            if not conn_str.endswith(".db"):
                raise ValueError("Invalid file type")
            self.conn = duckdb.connect(conn_str)
        elif conn_type.lower() == "postgres":
            pass
        elif conn_type.lower() == "mysql":
            pass
        else:
            raise ValueError("Invalid connection type")

    def get_full_table(self, table_name):
        return self.conn.execute(f"SELECT * FROM {table_name}").fetchdf()

    def from_csv(
        self,
        csv_path,
        table_name="temp_varDB",
    ):
        build_statement = (
            f"""CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM '{csv_path}'"""
        )
        logging.info(f"Building table {table_name} from {csv_path}")
        logging.debug(f"Executing: {build_statement}")
        self.conn.execute(build_statement)
        self.conn.execute(f'ALTER TABLE {table_name} DROP IF EXISTS index;')

    def to_csv(self, table_name, csv_path, split_column='CHROM'):
        if split_column:
            select_statement = f"SELECT DISTINCT({split_column}) FROM {table_name}"
            splits = self.conn.execute(select_statement).fetchall()
            for split in splits:
                logging.info(f"Writing file for {split[0]}")
                # Properly quote string values in the WHERE clause
                sql = f"""
                        SELECT * FROM {table_name} WHERE {split_column} = '{split[0]}'
                    """
                try:
                    self.conn.sql(sql).to_csv(f'{csv_path}_{split[0]}.csv')
                except Exception as e:
                    print(f"Error executing SQL: {e}")
        else:
            try:
                logging.info(f"Writing unified file for {table_name}")
                self.conn.execute(
                    f"COPY (SELECT * FROM {table_name}) TO '{csv_path}.csv'"
                )
            except Exception as e:
                print(f"Error executing SQL: {e}")

    def select_variants(
        self, table_name, chrom, start=None, end=None, ref=None, alt=None
    ):
        logging.info(f"Selecting variants for {chrom}")
        print(
            self.conn.execute(
                f"SELECT * FROM {table_name} WHERE CHROM = '{chrom}'"
            ).fetchone()
        )

    def select_variants_by_gene(self, table_name, gene):
        logging.info(f"Selecting variants for {gene}")
        print(
            self.conn.execute(
                f"SELECT * FROM {table_name} WHERE GENE = '{gene}'"
            ).fetchone()
        )

    def close(self):
        self.conn.close()


if __name__ == "__main__":
    vardb = VarDB()
    #    vardb.from_csv("/Users/juan.vazquez/dev/charly/data/variantsdb/*.csv")
    vardb.to_csv(
        "temp_varDB", "/Users/juan.vazquez/dev/charly/data/variantsdb/temp_varDB"
    )
