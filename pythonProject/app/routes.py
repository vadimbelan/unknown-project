# routes.py

import psycopg2
from flask import (Blueprint, request, redirect, url_for,
                   flash, session, render_template)
from functools import wraps
from werkzeug.security import generate_password_hash, check_password_hash
from .models import add_user, get_user_by_email
from datetime import datetime

from .scraper import parse_news
from .utils import get_location_info, write_html_for_two_cities, view_news_table, keywords
from .config import Config

sources = {
    1: {
        'name': 'Rambler.ru', 'table': 'news_nova_rambler_ru',
        'search_url': 'https://nova.rambler.ru/search?query={source}',
        'news_selector': 'h3.MixinCoolstream_news__title--1SXVJ', 'title_selector': None
    }
}

bp = Blueprint('main', __name__)

# Helper function to require login
def login_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if 'user_id' not in session:
            return redirect(url_for('main.login'))
        return f(*args, **kwargs)
    return decorated_function

@bp.route('/register', methods=['GET', 'POST'])
def register():
    if request.method == 'POST':
        username = request.form.get('username')
        email = request.form.get('email')
        password = request.form.get('password')
        hashed_password = generate_password_hash(password, method='pbkdf2:sha256')

        if get_user_by_email(email):
            flash("Email уже зарегистрирован!")
            return redirect(url_for('main.register'))

        # Добавляем пользователя в базу данных
        add_user(username, email, hashed_password)

        # Создаем сессию для пользователя
        user = get_user_by_email(email)
        session['user_id'] = user[0]
        session['username'] = user[1]

        # Перенаправляем на страницу задания бюджета после регистрации
        return redirect(url_for('main.budget'))

    return render_template('register.html')

@bp.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        email = request.form.get('email')
        password = request.form.get('password')

        user = get_user_by_email(email)
        if not user or not check_password_hash(user[3], password):
            flash("Неверный логин или пароль.")
            return redirect(url_for('main.login'))

        # Устанавливаем сессию для авторизованного пользователя
        session['user_id'] = user[0]
        session['username'] = user[1]
        flash("Успешный вход!")

        # Проверяем, есть ли у пользователя бюджет
        if not user[4]:  # Если бюджет отсутствует
            return redirect(url_for('main.budget'))

        return redirect(url_for('main.index'))

    return render_template('login.html')


@bp.route('/logout')
@login_required
def logout():
    session.pop('user_id', None)
    session.pop('username', None)
    flash("Вы вышли из системы.")
    return redirect(url_for('main.login'))

@bp.route('/')
@login_required
def index():
    user_id = session.get('user_id')

    conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")

    with conn:
        cur = conn.cursor()
        cur.execute("SELECT budget, currency FROM users WHERE id = %s", (user_id,))
        user_data = cur.fetchone()

    # Проверка, если данные пользователя не найдены
    if user_data is None:
        flash("Пользователь не найден.")
        return redirect(url_for('main.logout'))

    # Если бюджет не установлен, перенаправляем на страницу ввода бюджета
    if not user_data[0]:
        return redirect(url_for('main.budget'))

    budget, currency = user_data

    username = session.get('username')
    news_content = session.pop('news_content', '')

    return render_template('index.html', username=username, budget=budget, currency=currency, news_content=news_content)


@bp.route('/process', methods=['POST'])
def process():
    choice = int(request.form.get('source'))
    location_choice = request.form.get('location')
    use_location = location_choice == 'да'

    city = request.form.get('city')  # Получаем город, введенный пользователем

    tags = request.form.get('tags').strip().split('\n')
    tags = [tag.strip() for tag in tags if tag.strip()]

    if not city:
        flash("Необходимо указать город для отображения новостей.")
        return redirect(url_for('main.index'))

    selected_source = sources.get(choice, None)
    if not selected_source:
        flash("Неверный выбор. Пожалуйста, выберите источник из списка.")
        return redirect(url_for('main.index'))

    # Получаем информацию о введенном городе
    location_info_user_city = get_location_info_by_city(city)
    if not location_info_user_city:
        flash(f"Не удалось найти данные о городе {city}.")
        return redirect(url_for('main.index'))

    html_content = ""
    circle_color_user_city = '#000000'
    sentiment_message_user_city = "Нет данных о настроении."

    # Парсинг новостей для введенного города
    parse_news(
        city, location_info_user_city['Широта'], location_info_user_city['Долгота'],
        selected_source['table'], selected_source['search_url'],
        selected_source['news_selector'], selected_source['title_selector'], keywords
    )
    html_content += view_news_table(selected_source['table'], f"Новости из {city}")

    conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")

    sentiment_count_user_city = {'Positive': 0, 'Neutral': 0, 'Negative': 0}
    with conn:
        cur = conn.cursor()
        cur.execute(f"SELECT sentiment FROM {selected_source['table']}")
        sentiments_user_city= cur.fetchall()
    for sentiment in sentiments_user_city:
        sentiment_count_user_city[sentiment[0]] += 1


    max_sentiment_user_city = max(sentiment_count_user_city, key=sentiment_count_user_city.get, default='Neutral')
    color_map = {'Positive': '#00FF00', 'Neutral': '#0000FF', 'Negative': '#FF0000'}
    circle_color_user_city = color_map.get(max_sentiment_user_city, '#000000')
    sentiment_message_user_city = f"В {city} {max_sentiment_user_city.lower()} настроение."

    # Если выбран вывод по геопозиции, добавляем новости для текущего местоположения
    circle_color_geo = '#000000'
    sentiment_message_geo = "Нет данных о настроении."

    if use_location:
        location_info_geo = get_location_info(Config.DEFAULT_LOCATION)
        geo_city = location_info_geo['Город']

        # Парсинг новостей для геопозиции
        parse_news(
            geo_city, location_info_geo['Широта'], location_info_geo['Долгота'],
            selected_source['table'], selected_source['search_url'],
            selected_source['news_selector'], selected_source['title_selector'], keywords
        )
        html_content += view_news_table(selected_source['table'], f"Новости из {geo_city}")

        conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")

        sentiment_count_geo = {'Positive': 0, 'Neutral': 0, 'Negative': 0}
        with conn:
            cur = conn.cursor()
            cur.execute(f"SELECT sentiment FROM {selected_source['table']}")
            sentiments_geo = cur.fetchall()
        for sentiment in sentiments_geo:
            sentiment_count_geo[sentiment[0]] += 1


        max_sentiment_geo = max(sentiment_count_geo, key=sentiment_count_geo.get, default='Neutral')
        circle_color_geo = color_map.get(max_sentiment_geo, '#000000')
        sentiment_message_geo = f"В {geo_city} {max_sentiment_geo.lower()} настроение."

    # Создание карты с двумя точками
    html_output = write_html_for_two_cities(
        html_content, location_info_user_city, location_info_geo,
        selected_source['name'], datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        tags, circle_color_user_city, sentiment_message_user_city,
        circle_color_geo, sentiment_message_geo
    )

    session['news_content'] = html_output

    return redirect(url_for('main.index'))


def get_location_info_by_city(city):
    # Здесь можно использовать API для получения координат города (например, OpenCage или Google Geocoding)
    # Для примера, добавим несколько статических координат
    cities = {
        "Москва": {'Широта': 55.7558, 'Долгота': 37.6176, 'Страна': 'Россия', 'Регион': 'Москва', 'Город': 'Москва'},
        "Санкт-Петербург": {'Широта': 59.9343, 'Долгота': 30.3351, 'Страна': 'Россия', 'Регион': 'Санкт-Петербург', 'Город': 'Санкт-Петербург'}
    }
    return cities.get(city, None)


@bp.route('/add_favorite/<int:news_id>')
@login_required
def add_favorite(news_id):

    conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")

    with conn:
        cur = conn.cursor()
        cur.execute("INSERT INTO favorites (news_id) VALUES (%s)", (news_id,))
    return 'Success', 200

@bp.route('/remove_favorite/<int:news_id>')
@login_required
def remove_favorite(news_id):

    conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")

    with conn:
        cur = conn.cursor()
        cur.execute("DELETE FROM favorites WHERE news_id = %s", (news_id,))
    return 'Success', 200

@bp.route('/favorites')
@login_required
def favorites():

    conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")

    with conn:
        cur = conn.cursor()
        cur.execute('''SELECT n.id, n.title, n.latitude, n.longitude, n.sentiment, n.location
                       FROM news_nova_rambler_ru n
                       JOIN favorites f ON n.id = f.news_id''')
        favorite_news = cur.fetchall()

    return render_template('favorites.html', favorite_news=favorite_news)


@bp.route('/budget', methods=['GET', 'POST'])
@login_required
def budget():
    if request.method == 'POST':
        # Получаем данные бюджета и валюты
        budget = request.form.get('budget')
        currency = request.form.get('currency')

        user_id = session.get('user_id')

        conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")

        with conn:
            cur = conn.cursor()
            cur.execute("UPDATE users SET budget = %s, currency = %s WHERE id = %s", (budget, currency, user_id))

        # Перенаправляем пользователя на главную страницу
        return redirect(url_for('main.index'))

    return render_template('budget.html')
