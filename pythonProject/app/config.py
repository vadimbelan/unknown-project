#config.py

class Config:
    SECRET_KEY = 'your_secret_key'
    SQLALCHEMY_DATABASE_URI = 'postgresql://postgres:1234@localhost:5432/adventNews'
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    SESSION_TYPE = 'filesystem'
    SESSION_PERMANENT = False
    SESSION_USE_SIGNER = True
    SESSION_KEY_PREFIX = 'session:'
    SESSION_FILE_DIR = '/tmp/flask_session'
    SESSION_FILE_THRESHOLD = 100

    DEFAULT_LOCATION = {
        'Страна': 'Россия', 'Регион': 'Москва', 'Город': 'Нижний Новгород', 'Широта': 55.7558, 'Долгота': 37.6176
    }
