import os
import psycopg2

basedir = os.path.abspath(os.path.dirname(__file__))

def get_env_variable(name):
    try:
        return os.environ[name]
    except KeyError:
        message = "Expected environment variable '{}' not set.".format(name)
        raise Exception(message)

POSTGRES_URL="database-1.cmo9ouin7goc.us-east-1.rds.amazonaws.com"
POSTGRES_USER="postgres"
POSTGRES_PW="metabolomics2019"
POSTGRES_DB="postgres"

# the values of those depend on your setup
#POSTGRES_URL = get_env_variable("POSTGRES_URL")
#POSTGRES_USER = get_env_variable("POSTGRES_USER")
#POSTGRES_PW = get_env_variable("POSTGRES_PW")
#POSTGRES_DB = get_env_variable("POSTGRES_DB")

DB_URL = 'postgresql+psycopg2://{user}:{pw}@{url}/{db}'.format(user=POSTGRES_USER,pw=POSTGRES_PW,url=POSTGRES_URL,db=POSTGRES_DB)

#'DATABASE_URL'
class BaseConfig:
    SQLALCHEMY_DATABASE_URI = (DB_URL)
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    SECRET_KEY = 'secret_key'
    #SECRET_KEY = os.environ['SECRET_KEY']
