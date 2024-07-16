import requests
from bs4 import BeautifulSoup
import geocoder
import pymorphy2
import sqlite3
from googletrans import Translator
from transformers import pipeline
from datetime import datetime

default_location = {
    'Страна': '', 'Регион': '', 'Город': 'Нижний Новгород', 'Широта': 0.0, 'Долгота': 0.0
}

translations = {'country': {}, 'region': {}, 'city': {}}
keywords = {'positive': ["Площадки для кинопоказов под открытым небом заработали"], 'neutral': [], 'negative': []}

morph = pymorphy2.MorphAnalyzer()
translator = Translator()
sentiment_analysis = pipeline("sentiment-analysis", model="distilbert-base-uncased-finetuned-sst-2-english")

def init_db():
    with sqlite3.connect('news.db') as conn:
        cursor = conn.cursor()
        for table in ['news_mail_ru', 'news_ria_ru', 'news_nova_rambler_ru']:
            cursor.execute(f"DROP TABLE IF EXISTS {table}")
            cursor.execute(f'''CREATE TABLE {table} (
                                id INTEGER PRIMARY KEY AUTOINCREMENT,
                                title TEXT,
                                latitude REAL,
                                longitude REAL,
                                sentiment TEXT)''')

def add_news_to_db(table, title, latitude, longitude, sentiment):
    with sqlite3.connect('news.db') as conn:
        cursor = conn.cursor()
        cursor.execute(f"INSERT INTO {table} (title, latitude, longitude, sentiment) VALUES (?, ?, ?, ?)",
                       (title, latitude, longitude, sentiment))
        cursor.execute(f"DELETE FROM {table} WHERE id NOT IN (SELECT id FROM {table} ORDER BY id DESC LIMIT 20)")

def get_location_info(defaults):
    try:
        g = geocoder.ip('me')
        if g.ok:
            return {
                'Страна': defaults['Страна'] if defaults['Страна'] else translations['country'].get(g.country, g.country),
                'Регион': defaults['Регион'] if defaults['Регион'] else translations['region'].get(g.state, g.state),
                'Город': defaults['Город'] if defaults['Город'] else translations['city'].get(g.city, g.city),
                'Широта': defaults['Широта'] if defaults['Широта'] else g.latlng[0],
                'Долгота': defaults['Долгота'] if defaults['Долгота'] else g.latlng[1]
            }
    except Exception as e:
        defaults['Ошибка'] = str(e)
    return defaults

def translate_text(text, dest_language='en'):
    try:
        return translator.translate(text, dest=dest_language).text
    except Exception as e:
        return f"Translation error: {e}"

def get_sentiment(text):
    for sentiment, kw_list in keywords.items():
        if any(kw in text for kw in kw_list):
            return sentiment.capitalize()

    sentiment_result = sentiment_analysis(translate_text(text))[0]
    return 'Neutral' if sentiment_result['score'] < 0.9 else sentiment_result['label'].capitalize()

def parse_news(source, latitude, longitude, table, search_url, news_selector, title_selector):
    response = requests.get(search_url.format(source=source))
    soup = BeautifulSoup(response.content, 'html.parser')

    news_list = []
    for item in soup.select(news_selector)[:10]:
        title = item.select_one(title_selector).get_text(strip=True) if title_selector else item.get_text(strip=True)
        sentiment = get_sentiment(title)
        add_news_to_db(table, title, latitude, longitude, sentiment)
        news_list.append({'title': title, 'latitude': latitude, 'longitude': longitude, 'sentiment': sentiment})
        print(f"Original: {title}\nTranslated: {translate_text(title)}\nSentiment: {sentiment}\n")
    return news_list

def view_news_table(table_name):
    with sqlite3.connect('news.db') as conn:
        cursor = conn.cursor()
        rows = cursor.execute(f"SELECT * FROM {table_name}").fetchall()

    table_html = f"<h2>{table_name.replace('_', ' ').title()} Table</h2>\n<table border='1'>\n"
    table_html += "<tr><th>ID</th><th>Title</th><th>Latitude</th><th>Longitude</th><th>Sentiment</th></tr>\n"
    for row in rows:
        color = {'Positive': 'green', 'Neutral': 'blue', 'Negative': 'red'}.get(row[4], 'black')
        table_html += f"<tr><td>{row[0]}</td><td>{row[1]}</td><td>{row[2]}</td><td>{row[3]}</td><td style='color:{color}'>{row[4]}</td></tr>\n"
    return table_html + "</table>\n"

def write_html(content, location_info, source, timestamp):
    with open("output.html", "w", encoding="utf-8") as html_file:
        html_file.write(f"<html><head><title>News Output</title><meta charset='UTF-8'></head><body>")
        html_file.write(f"<h2>Source Information</h2><p>Источник: {source}</p>")
        html_file.write(f"<p>Местоположение: {location_info}</p>")
        html_file.write(f"<p>Время: {timestamp}</p>{content}</body></html>")

timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
html_content = ""

location_info = get_location_info(default_location)
source = location_info.get('Город', '')

sources = {
    1: {'name': 'Rambler.ru', 'table': 'news_nova_rambler_ru', 'search_url': 'https://nova.rambler.ru/search?query={source}', 'news_selector': 'h3.MixinCoolstream_news__title--1SXVJ', 'title_selector': None}
}

init_db()

print("Список новостей:")
for number, info in sources.items():
    print(f"{number}. {info['name']}")

try:
    choice = int(input("Введите номер: "))
    if choice in sources:
        src = sources[choice]
        parse_news(source, location_info['Широта'], location_info['Долгота'], src['table'], src['search_url'], src['news_selector'], src['title_selector'])
        html_content += view_news_table(src['table'])
    else:
        print("Неверный выбор. Пожалуйста, выберите номер от 1 до 3.")
except ValueError:
    print("Пожалуйста, введите номер из списка.")

write_html(html_content, location_info, source, timestamp)
