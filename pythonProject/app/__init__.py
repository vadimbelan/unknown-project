# __init__.py

from flask import Flask
from flask_session import Session  # Import Flask-Session

def create_app():
    app = Flask(__name__)
    app.config.from_object('app.config.Config')

    # Initialize Flask-Session
    Session(app)

    from . import routes
    app.register_blueprint(routes.bp)

    return app
