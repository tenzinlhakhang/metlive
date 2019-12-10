
import psycopg2



engine = psycopg2.connect(
    database="postgres",
    user="Tenzin",
    password="tenzin2014",
    host="database-3.cmo9ouin7goc.us-east-1.rds.amazonaws.com",
    port='5432'
)


