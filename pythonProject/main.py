import requests
from bs4 import BeautifulSoup
import geocoder

default_location = {
    'Страна': '',
    'Регион': '',
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
def parse_news(source):
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
        news_list.append({'title': title, 'link': link, 'date': date})

    return news_list


# Сохраняем информацию о местоположении в переменную location_info
location_info = get_location_info(default_location)

# Выводим информацию о местоположении на экран
print(location_info)

# Фильтрация местоположения по конкретной области
source = location_info.get('Регион', '')

# Парсим новости
if source:
    news = parse_news(source)
    for news_item in news:
        print(f"Title: {news_item['title']}\nLink: {news_item['link']}\nDate: {news_item['date']}\n")
else:
    print("Не удалось определить местоположение")
