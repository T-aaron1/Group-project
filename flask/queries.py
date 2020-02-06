# close db

import sqlite3
import pandas as pd

def select_gral(DATABASE, select_text,from_table, where_data):
    query = "SELECT {} FROM {} WHERE {}".format(select_text,from_table, where_data)
    db = sqlite3.connect(DATABASE)
    db_results = pd.read_sql_query(query, db)
    db.close()
    return  db_results


def query_n_results(DATABASE, select_text,from_table, where_data):
    query = "SELECT COUNT({}) FROM {} WHERE {}".format(select_text,from_table, where_data)
    db = sqlite3.connect(DATABASE)
    n_results = pd.read_sql_query(query, db).iloc[0,0]
    db.close()
    return  n_results


def query_is_unique(DATABASE, select_text,from_table, where_data):
    n_results = query_n_results(DATABASE, select_text,from_table, where_data)
    return  n_results == 1
