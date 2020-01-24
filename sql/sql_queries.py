import sqlite3
import pandas as pd

db = sqlite3.connect('/homes/dtg30/Desktop/group_proj_2/kinase_project.db')
cursor = db.cursor()

cursor.execute("SELECT * from families")
table = pd.read_sql_query("SELECT * from families", db)
