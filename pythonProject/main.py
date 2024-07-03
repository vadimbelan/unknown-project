import requests
from bs4 import BeautifulSoup
import geocoder
import pymorphy2
import sqlite3

# Инициализация значений по умолчанию для местоположения
default_location = {
    'Страна': '',
    'Регион': 'Печора',
    'Город': '',
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


# Создание таблицы новостей
def create_and_clear_news_table():
    conn = sqlite3.connect('news.db')
    cursor = conn.cursor()
    cursor.execute('''DROP TABLE IF EXISTS news''')  # Удаление существующей таблицы
    cursor.execute('''CREATE TABLE news (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        title TEXT,
                        link TEXT,
                        date TEXT,
                        latitude REAL,
                        longitude REAL
                      )''')
    conn.commit()
    conn.close()


create_and_clear_news_table()


# Функция для добавления новостей в базу данных
def add_news_to_db(title, link, date, latitude, longitude):
    conn = sqlite3.connect('news.db')
    cursor = conn.cursor()
    cursor.execute("INSERT INTO news (title, link, date, latitude, longitude) VALUES (?, ?, ?, ?, ?)",
                   (title, link, date, latitude, longitude))
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

        # Используем координаты из get_location_info
        add_news_to_db(title, link, date, latitude, longitude)
        news_list.append({'title': title, 'link': link, 'date': date, 'latitude': latitude, 'longitude': longitude})

    return news_list


# Функция для лемматизации текста
def lemmatize_text(text):
    words = text.split()
    lemmas = [morph.parse(word)[0].normal_form for word in words]
    return ' '.join(lemmas)


# Функция для вывода реакции на новости
def react_to_news(news_title):
    reactions = {
        'Ключ': 'Значение',
    }

    # Лемматизация заголовка новости
    lemmatized_title = lemmatize_text(news_title.lower())

    for keyword, reaction in reactions.items():
        if keyword in lemmatized_title:
            return reaction
    return 'Прочитайте последние новости!'


# Функция для отображения таблицы новостей
def view_news_table():
    conn = sqlite3.connect('news.db')
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM news")
    rows = cursor.fetchall()
    conn.close()

    for row in rows:
        print(row)


# Сохраняем информацию о местоположении в переменную location_info
location_info = get_location_info(default_location)

# Выводим информацию о местоположении на экран
print(location_info)

# Фильтрация местоположения по конкретной области
source = location_info.get('Город', '')

# Парсим новости
if source:
    news = parse_news(source, location_info['Широта'], location_info['Долгота'])
    for news_item in news:
        reaction = react_to_news(news_item['title'])
        print(
            f"Title: {news_item['title']}\nLink: {news_item['link']}\nDate: {news_item['date']}\nReaction: {reaction}\n")
else:
    print("Не удалось определить местоположение")

# Вызов функции для отображения таблицы
view_news_table()
