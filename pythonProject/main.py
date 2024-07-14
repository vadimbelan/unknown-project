import requests
from bs4 import BeautifulSoup
import geocoder
import pymorphy2
import sqlite3
from googletrans import Translator
from transformers import pipeline
from datetime import datetime

default_location = {
    'Страна': '',
    'Регион': '',
    'Город': 'Сочи',
    'Широта': 0.0,
    'Долгота': 0.0
}

country_translations = {
    '': '',
}

region_translations = {
    '': '',
}

city_translations = {
    '': '',
}

positive_keywords = []
neutral_keywords = ["Главные темы часа"]
negative_keywords = []

morph = pymorphy2.MorphAnalyzer()
translator = Translator()
sentiment_analysis = pipeline("sentiment-analysis", model="distilbert-base-uncased-finetuned-sst-2-english")


def create_and_clear_news_tables():
    conn = sqlite3.connect('news.db')
    cursor = conn.cursor()
    cursor.execute('''DROP TABLE IF EXISTS news_mail_ru''')
    cursor.execute('''DROP TABLE IF EXISTS news_ria_ru''')
    cursor.execute('''DROP TABLE IF EXISTS news_nova_rambler_ru''')
    cursor.execute('''CREATE TABLE news_mail_ru (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        title TEXT,
                        link TEXT,
                        latitude REAL,
                        longitude REAL,
                        sentiment TEXT
                      )''')
    cursor.execute('''CREATE TABLE news_ria_ru (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        title TEXT,
                        link TEXT,
                        date TEXT,
                        latitude REAL,
                        longitude REAL,
                        sentiment TEXT
                      )''')
    cursor.execute('''CREATE TABLE news_nova_rambler_ru (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        title TEXT,
                        latitude REAL,
                        longitude REAL,
                        sentiment TEXT
                      )''')
    conn.commit()
    conn.close()


create_and_clear_news_tables()


def add_news_to_db(table, title, link, date, latitude, longitude, sentiment):
    conn = sqlite3.connect('news.db')
    cursor = conn.cursor()
    if table == 'news_mail_ru':
        cursor.execute(f"INSERT INTO {table} (title, link, latitude, longitude, sentiment) VALUES (?, ?, ?, ?, ?)",
                       (title, link, latitude, longitude, sentiment))
    elif table == 'news_nova_rambler_ru':
        cursor.execute(f"INSERT INTO {table} (title, latitude, longitude, sentiment) VALUES (?, ?, ?, ?)",
                       (title, latitude, longitude, sentiment))
    else:
        cursor.execute(
            f"INSERT INTO {table} (title, link, date, latitude, longitude, sentiment) VALUES (?, ?, ?, ?, ?, ?)",
            (title, link, date, latitude, longitude, sentiment))

    # Удаление лишних записей для соблюдения лимита в 5 строк
    cursor.execute(f"DELETE FROM {table} WHERE id NOT IN (SELECT id FROM {table} ORDER BY id DESC LIMIT 5)")
    conn.commit()
    conn.close()


def get_location_info(defaults):
    try:
        g = geocoder.ip('me')

        if g.ok:
            country = g.country if defaults['Страна'] == '' else defaults['Страна']
            region = g.state if defaults['Регион'] == '' else defaults['Регион']
            city = g.city if defaults['Город'] == '' else defaults['Город']

            return {
                'Страна': country_translations.get(country, country),
                'Регион': region_translations.get(region, region),
                'Город': city_translations.get(city, city),
                'Широта': g.latlng[0] if defaults['Широта'] == 0.0 else defaults['Широта'],
                'Долгота': g.latlng[1] if defaults['Долгота'] == 0.0 else defaults['Долгота']
            }
        else:
            return defaults
    except Exception as e:
        defaults['Ошибка'] = str(e)
        return defaults


def translate_text(text, dest_language='en'):
    try:
        translated = translator.translate(text, dest=dest_language)
        return translated.text
    except Exception as e:
        return f"Translation error: {e}"


def get_sentiment(text):
    for keyword in positive_keywords:
        if keyword in text:
            return 'Positive'
    for keyword in neutral_keywords:
        if keyword in text:
            return 'Neutral'
    for keyword in negative_keywords:
        if keyword in text:
            return 'Negative'

    translated_text = translate_text(text)
    sentiment_result = sentiment_analysis(translated_text)[0]
    sentiment_score = sentiment_result['score']
    if sentiment_score < 0.9:
        return 'Neutral'
    elif sentiment_result['label'] == 'POSITIVE':
        return 'Positive'
    else:
        return 'Negative'


def parse_news_mail_ru(source, latitude, longitude):
    search_url = f'https://news.mail.ru/search/?q={source}'
    response = requests.get(search_url)
    soup = BeautifulSoup(response.content, 'html.parser')

    news_items = soup.find_all('a', class_='newsitem__title link-holder')

    news_list = []
    for item in news_items[:5]:  # Ограничение количества обрабатываемых новостей до 5
        title = item.get_text(strip=True)
        link = item['href']

        sentiment = get_sentiment(title)

        add_news_to_db('news_mail_ru', title, link, None, latitude, longitude, sentiment)
        news_list.append(
            {'title': title, 'link': link, 'latitude': latitude, 'longitude': longitude, 'sentiment': sentiment})

        translated_title = translate_text(title)
        print(f"Original: {title}\nTranslated: {translated_title}\nSentiment: {sentiment}\n")

    return news_list


def parse_news_ria_ru(source, latitude, longitude):
    search_url = f'https://ria.ru/search/?query={source}'
    response = requests.get(search_url)
    soup = BeautifulSoup(response.content, 'html.parser')

    news_items = soup.find_all('div', class_='list-item')

    news_list = []
    for item in news_items[:5]:  # Ограничение количества обрабатываемых новостей до 5
        title = item.find('a', class_='list-item__title').get_text(strip=True)
        link = item.find('a', class_='list-item__title')['href']
        date = item.find('div', class_='list-item__date').get_text(strip=True)

        sentiment = get_sentiment(title)

        add_news_to_db('news_ria_ru', title, link, date, latitude, longitude, sentiment)
        news_list.append({'title': title, 'link': link, 'date': date, 'latitude': latitude, 'longitude': longitude,
                          'sentiment': sentiment})

        translated_title = translate_text(title)
        print(f"Original: {title}\nTranslated: {translated_title}\nSentiment: {sentiment}\n")

    return news_list


def parse_news_nova_rambler_ru(source, latitude, longitude):
    search_url = f'https://nova.rambler.ru/search?query={source}'
    response = requests.get(search_url)
    soup = BeautifulSoup(response.content, 'html.parser')

    news_items = soup.find_all('h3', class_="MixinCoolstream_news__title--1SXVJ")

    news_list = []
    for item in news_items[:5]:  # Ограничение количества обрабатываемых новостей до 5
        title = item.get_text(strip=True)

        sentiment = get_sentiment(title)

        add_news_to_db('news_nova_rambler_ru', title, None, None, latitude, longitude, sentiment)
        news_list.append({'title': title, 'latitude': latitude, 'longitude': longitude, 'sentiment': sentiment})

        translated_title = translate_text(title)
        print(f"Original: {title}\nTranslated: {translated_title}\nSentiment: {sentiment}\n")

    return news_list


def view_news_table(table_name):
    conn = sqlite3.connect('news.db')
    cursor = conn.cursor()
    cursor.execute(f"SELECT * FROM {table_name}")
    rows = cursor.fetchall()
    conn.close()

    table_html = f"<h2>{table_name.replace('_', ' ').title()} Table</h2>\n<table border='1'>\n"
    table_html += "<tr><th>ID</th><th>Title</th><th>Link</th>"
    if table_name == 'news_ria_ru':
        table_html += "<th>Date</th>"
    table_html += "<th>Latitude</th><th>Longitude</th><th>Sentiment</th></tr>\n"
    for row in rows:
        sentiment_color = {
            'Positive': 'green',
            'Neutral': 'blue',
            'Negative': 'red'
        }.get(row[5 if table_name == 'news_mail_ru' else 6 if table_name == 'news_ria_ru' else 4], 'black')
        table_html += f"<tr><td>{row[0]}</td><td>{row[1]}</td><td>{row[2] if table_name != 'news_nova_rambler_ru' else 'N/A'}</td>"
        if table_name == 'news_ria_ru':
            table_html += f"<td>{row[3]}</td>"
        table_html += f"<td>{row[3 if table_name == 'news_mail_ru' else 4 if table_name == 'news_ria_ru' else 2]}</td><td>{row[4 if table_name == 'news_mail_ru' else 5 if table_name == 'news_ria_ru' else 3]}</td><td style='color:{sentiment_color}'>{row[5 if table_name == 'news_mail_ru' else 6 if table_name == 'news_ria_ru' else 4]}</td></tr>\n"
    table_html += "</table>\n"

    return table_html


def write_html(content, location_info, source, timestamp):
    with open("output.html", "w", encoding="utf-8") as html_file:
        html_file.write("<html>\n<head>\n<title>News Output</title>\n<meta charset='UTF-8'>\n</head>\n<body>\n")
        html_file.write("<h1>News Output</h1>\n")
        html_file.write(f"<h2>Source Information</h2>\n<p>Новости взяты с источника: {source}</p>\n")
        html_file.write(f"<p>Новости парсятся по местоположению: {location_info}</p>\n")
        html_file.write(f"<p>Время парсинга: {timestamp}</p>\n")
        html_file.write(content)
        html_file.write("</body>\n</html>")


# Get current timestamp
timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

html_content = ""

location_info = get_location_info(default_location)

source = location_info.get('Город', '')

sources = {
    1: 'Mail.ru',
    2: 'RIA.ru',
    3: 'Rambler.ru'
}

print("Выберите источник новостей:")
for number, name in sources.items():
    print(f"{number}. {name}")

try:
    choice = int(input("Введите номер источника: "))
    if choice in sources:
        if choice == 1:
            parse_news_mail_ru(source, location_info['Широта'], location_info['Долгота'])
            html_content += view_news_table('news_mail_ru')
        elif choice == 2:
            parse_news_ria_ru(source, location_info['Широта'], location_info['Долгота'])
            html_content += view_news_table('news_ria_ru')
        elif choice == 3:
            parse_news_nova_rambler_ru(source, location_info['Широта'], location_info['Долгота'])
            html_content += view_news_table('news_nova_rambler_ru')
    else:
        print("Неверный выбор. Пожалуйста, выберите номер от 1 до 3.")
except ValueError:
    print("Пожалуйста, введите номер источника.")

# Write HTML with detailed information
write_html(html_content, location_info, source, timestamp)
