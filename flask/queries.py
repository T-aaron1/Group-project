# close db

import sqlite3
import pandas as pd

def select_gral(DATABASE, select_text,from_table, where_data, like_something):
    query = "SELECT {} FROM {} WHERE {} LIKE '{}'".format(select_text,from_table, where_data, like_something)
    db = sqlite3.connect(DATABASE)
    db_results = pd.read_sql_query(query, db)
    db.close()
    return  db_results

def query_n_results(DATABASE, select_text,from_table, where_data, like_something):
    db_results = select_gral(DATABASE, select_text,from_table, where_data, like_something)
    n_results = db_results.shape[0]
    return  n_results


def query_is_unique(DATABASE, select_text,from_table, where_data, like_something):
    n_results = query_n_results(DATABASE, select_text,from_table, where_data, like_something)
    return  n_results == 1
