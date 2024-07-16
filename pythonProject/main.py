import requests
from bs4 import BeautifulSoup
import geocoder
import pymorphy2
import sqlite3
from googletrans import Translator
from transformers import pipeline
from datetime import datetime
import spacy

default_location = {
    'Страна': '', 'Регион': '', 'Город': 'Нижний Новгород', 'Широта': 0.0, 'Долгота': 0.0
}

translations = {
    'country': {}, 'region': {}, 'city': {}
}

keywords = {
    'positive': [],
    'neutral': [],
    'negative': []
}

morph = pymorphy2.MorphAnalyzer()
translator = Translator()
model1 = pipeline("sentiment-analysis", model="nlptown/bert-base-multilingual-uncased-sentiment")
model2 = pipeline("sentiment-analysis", model="distilbert-base-uncased-finetuned-sst-2-english")
nlp = spacy.load("ru_core_news_sm")


def init_db():
    with sqlite3.connect('news.db') as conn:
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS news_nova_rambler_ru")
        cur.execute('''CREATE TABLE news_nova_rambler_ru (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        title TEXT,
                        latitude REAL,
                        longitude REAL,
                        sentiment TEXT,
                        location TEXT)''')
        cur.execute("DROP TABLE IF EXISTS news_tags")
        cur.execute('''CREATE TABLE news_tags (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        title TEXT,
                        latitude REAL,
                        longitude REAL,
                        sentiment TEXT,
                        location TEXT)''')


def add_news_to_db(table, title, latitude, longitude, sentiment, location):
    with sqlite3.connect('news.db') as conn:
        cur = conn.cursor()
        cur.execute(f"INSERT INTO {table} (title, latitude, longitude, sentiment, location) VALUES (?, ?, ?, ?, ?)",
                    (title, latitude, longitude, sentiment, location))
        cur.execute(f"DELETE FROM {table} WHERE id NOT IN (SELECT id FROM {table} ORDER BY id DESC LIMIT 20)")


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

    sentiment_result1 = model1(text)[0]
    translated_text = translate_text(text)
    sentiment_result2 = model2(translated_text)[0]

    label_map1 = {
        '1 star': 'Negative', '2 stars': 'Negative', '3 stars': 'Neutral',
        '4 stars': 'Positive', '5 stars': 'Positive'
    }
    sentiment1 = label_map1.get(sentiment_result1['label'], 'Neutral')

    sentiment2 = 'Negative' if sentiment_result2['label'] == 'NEGATIVE' else 'Positive' if sentiment_result2['label'] == 'POSITIVE' else 'Neutral'

    print(f"Original text: {text}")
    print(f"Translated text: {translated_text}")
    print(f"Sentiment by model 1 (Russian): {sentiment1}")
    print(f"Sentiment by model 2 (English): {sentiment2}")

    final_sentiment = sentiment1 if sentiment1 == sentiment2 else 'Neutral'

    print(f"Final sentiment: {final_sentiment}\n")
    return final_sentiment


def extract_location(text):
    doc = nlp(text)
    locations = [ent.text for ent in doc.ents if ent.label_ == "LOC"]
    if locations:
        return ', '.join(locations)
    return None


def parse_news(source, latitude, longitude, table, search_url, news_selector, title_selector):
    response = requests.get(search_url.format(source=source))
    soup = BeautifulSoup(response.content, 'html.parser')

    news_list = []
    for item in soup.select(news_selector)[:10]:
        title = item.select_one(title_selector).get_text(strip=True) if title_selector else item.get_text(strip=True)
        sentiment = get_sentiment(title)
        location = extract_location(title)
        add_news_to_db(table, title, latitude, longitude, sentiment, location)
        news_list.append({
            'title': title, 'latitude': latitude, 'longitude': longitude, 'sentiment': sentiment, 'location': location
        })
    return news_list


def view_news_table(table_name, title):
    with sqlite3.connect('news.db') as conn:
        cur = conn.cursor()
        rows = cur.execute(f"SELECT * FROM {table_name}").fetchall()

    table_html = f"<h2>{title}</h2>\n<table border='1'>\n"
    table_html += "<tr><th>ID</th><th>Title</th><th>Latitude</th><th>Longitude</th><th>Sentiment</th><th>Location</th></tr>\n"
    for row in rows:
        color = {'Positive': 'green', 'Neutral': 'blue', 'Negative': 'red'}.get(row[4], 'black')
        table_html += f"<tr><td>{row[0]}</td><td>{row[1]}</td><td>{row[2]}</td><td>{row[3]}</td><td style='color:{color}'>{row[4]}</td><td>{row[5]}</td></tr>\n"
    return table_html + "</table>\n"


def write_html(content, loc_info, src, ts, tags):
    with open("output.html", "w", encoding="utf-8") as html_file:
        html_file.write(f"<html><head><title>News Output</title><meta charset='UTF-8'></head><body>")
        html_file.write(f"<h2>Source Information</h2><p>Источник: {src}</p>")
        html_file.write(f"<p>Местоположение: {loc_info}</p>")
        html_file.write(f"<p>Время: {ts}</p>")
        html_file.write(f"<h2>Введённые теги:</h2><p>{', '.join(tags)}</p>")
        html_file.write(content + "</body></html>")


timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
html_content = ""

location_info = get_location_info(default_location)

sources = {
    1: {
        'name': 'Rambler.ru', 'table': 'news_nova_rambler_ru',
        'search_url': 'https://nova.rambler.ru/search?query={source}',
        'news_selector': 'h3.MixinCoolstream_news__title--1SXVJ', 'title_selector': None
    }
}

init_db()

print("Список новостей:")
for number, info in sources.items():
    print(f"{number}. {info['name']}")

try:
    choice = int(input("Введите номер: "))
    if choice not in sources:
        print("Неверный выбор. Пожалуйста, выберите номер из списка.")
        exit()

    location_choice = input("Выводить новости по геопозиции? (да/нет): ").strip().lower()
    use_location = location_choice == 'да'

    tags = []
    print("Введите теги (каждый тег - через enter, чтобы закончить вводить теги - введите 0): ", end='')
    while True:
        tag = input().strip()
        if tag == '0':
            break
        if tag:
            tags.append(tag)
        print(f"Введённые теги: {', '.join(tags)}")
        print("Введите теги (каждый тег - через enter, чтобы закончить вводить теги - введите 0): ", end='')

    if not use_location and not tags:
        print("Не удалось вывести последние новости, недостаточно входных данных.")
        exit()

    selected_source = sources[choice]

    if use_location:
        chosen_source = location_info.get('Город', '')
        parse_news(
            chosen_source, location_info['Широта'], location_info['Долгота'],
            selected_source['table'], selected_source['search_url'],
            selected_source['news_selector'], selected_source['title_selector']
        )
        html_content += view_news_table(selected_source['table'], "News by Geolocation")

    if tags:
        with sqlite3.connect('news.db') as conn:
            cur = conn.cursor()
            cur.execute("DELETE FROM news_tags")

        for tag in tags:
            parse_news(
                tag, location_info['Широта'], location_info['Долгота'], 'news_tags',
                selected_source['search_url'], selected_source['news_selector'], selected_source['title_selector']
            )
        html_content += view_news_table('news_tags', "News by Tags")

    write_html(html_content, location_info, selected_source['name'], timestamp, tags)

except ValueError:
    print("Пожалуйста, введите номер из списка.")
