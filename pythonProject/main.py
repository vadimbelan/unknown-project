import requests
from bs4 import BeautifulSoup
import geocoder
import pymorphy2
import sqlite3
from googletrans import Translator
from transformers import pipeline

# Инициализация значений по умолчанию для местоположения
default_location = {
    'Страна': '',
    'Регион': '',
    'Город': 'Печора',
    'Широта': 0.0,
    'Долгота': 0.0
}

# Словари для перевода стран, регионов и городов
country_translations = {
    'RU': 'Россия',
}

region_translations = {
    'Komi': 'Коми',
}

city_translations = {
    'Yemva': 'Емва',
}

# Инициализация анализатора pymorphy2
morph = pymorphy2.MorphAnalyzer()

# Инициализация переводчика
translator = Translator()

# Инициализация анализатора тональности на базе DistilBERT
sentiment_analysis = pipeline("sentiment-analysis", model="distilbert-base-uncased-finetuned-sst-2-english")

# Создание таблиц новостей
def create_and_clear_news_tables():
    conn = sqlite3.connect('news.db')
    cursor = conn.cursor()
    cursor.execute('''DROP TABLE IF EXISTS news''')
    cursor.execute('''CREATE TABLE news (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        title TEXT,
                        link TEXT,
                        date TEXT,
                        latitude REAL,
                        longitude REAL,
                        sentiment TEXT
                      )''')
    cursor.execute('''DROP TABLE IF EXISTS positive_news''')
    cursor.execute('''CREATE TABLE positive_news AS SELECT * FROM news WHERE 0''')
    cursor.execute('''DROP TABLE IF EXISTS negative_news''')
    cursor.execute('''CREATE TABLE negative_news AS SELECT * FROM news WHERE 0''')
    cursor.execute('''DROP TABLE IF EXISTS neutral_news''')
    cursor.execute('''CREATE TABLE neutral_news AS SELECT * FROM news WHERE 0''')
    conn.commit()
    conn.close()

create_and_clear_news_tables()

# Функция для добавления новостей в базу данных
def add_news_to_db(title, link, date, latitude, longitude, sentiment):
    conn = sqlite3.connect('news.db')
    cursor = conn.cursor()
    cursor.execute("INSERT INTO news (title, link, date, latitude, longitude, sentiment) VALUES (?, ?, ?, ?, ?, ?)",
                   (title, link, date, latitude, longitude, sentiment))
    conn.commit()
    cursor.execute(f"INSERT INTO {sentiment.lower()}_news SELECT * FROM news WHERE id = (SELECT MAX(id) FROM news)")
    conn.commit()
    conn.close()

# Получаем информацию о местоположении
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

# Функция для перевода текста на английский
def translate_text(text, dest_language='en'):
    try:
        translated = translator.translate(text, dest=dest_language)
        return translated.text
    except Exception as e:
        return f"Translation error: {e}"

# Функция для определения тональности текста
def get_sentiment(text):
    translated_text = translate_text(text)
    sentiment_result = sentiment_analysis(translated_text)[0]
    sentiment_score = sentiment_result['score']
    if sentiment_score < 0.7:
        return 'Neutral'
    elif sentiment_result['label'] == 'POSITIVE':
        return 'Positive'
    else:
        return 'Negative'

# Функция для парсинга новостей
def parse_news(source, latitude, longitude):
    search_url = f'https://ria.ru/search/?query={source}'
    response = requests.get(search_url)
    soup = BeautifulSoup(response.content, 'html.parser')

    # Найти все элементы, содержащие новости
    news_items = soup.find_all('div', class_='list-item')

    news_list = []
    for item in news_items:
        title = item.find('a', class_='list-item__title').get_text(strip=True)
        link = item.find('a', class_='list-item__title')['href']
        date = item.find('div', class_='list-item__date').get_text(strip=True)

        sentiment = get_sentiment(title)

        # Используем координаты из get_location_info
        add_news_to_db(title, link, date, latitude, longitude, sentiment)
        news_list.append({'title': title, 'link': link, 'date': date, 'latitude': latitude, 'longitude': longitude, 'sentiment': sentiment})

        # Перевод заголовка и вывод в консоль
        translated_title = translate_text(title)
        print(f"Original: {title}\nTranslated: {translated_title}\nSentiment: {sentiment}\n")

    return news_list

# Функция для отображения таблицы новостей
def view_news_table(table_name):
    conn = sqlite3.connect('news.db')
    cursor = conn.cursor()
    cursor.execute(f"SELECT * FROM {table_name}")
    rows = cursor.fetchall()
    conn.close()

    table_html = f"<h2>{table_name.replace('_', ' ').title()} Table</h2>\n<table border='1'>\n"
    table_html += "<tr><th>ID</th><th>Title</th><th>Link</th><th>Date</th><th>Latitude</th><th>Longitude</th><th>Sentiment</th></tr>\n"
    for row in rows:
        table_html += f"<tr><td>{row[0]}</td><td>{row[1]}</td><td><a href='{row[2]}'>{row[2]}</a></td><td>{row[3]}</td><td>{row[4]}</td><td>{row[5]}</td><td>{row[6]}</td></tr>\n"
    table_html += "</table>\n"

    return table_html

# Функция для записи HTML-файла
def write_html(content):
    with open("output.html", "w", encoding="utf-8") as html_file:
        html_file.write("<html>\n<head>\n<title>News Output</title>\n<meta charset='UTF-8'>\n</head>\n<body>\n")
        html_file.write(content)
        html_file.write("</body>\n</html>")

# Собираем весь HTML-контент
html_content = ""

# Сохраняем информацию о местоположении в переменную location_info
location_info = get_location_info(default_location)

# Добавляем информацию о местоположении в HTML-контент
html_content += "<h1>Location Info</h1>\n"
html_content += f"<p>{location_info}</p>\n"

# Фильтрация местоположения по конкретной области
source = location_info.get('Город', '')

# Парсим новости
if source:
    parse_news(source, location_info['Широта'], location_info['Долгота'])
else:
    html_content += "<p>Не удалось определить местоположение</p>\n"

# Добавляем таблицы новостей в HTML-контент
html_content += view_news_table('news')
html_content += view_news_table('positive_news')
html_content += view_news_table('negative_news')
html_content += view_news_table('neutral_news')

# Запись финального HTML-файла
write_html(html_content)
