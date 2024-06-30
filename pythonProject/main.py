import geocoder

# Заглушки для местоположения
default_location = {
    'Страна': '',
    'Регион': '',
    'Город': '',
    'Широта': 0.0,
    'Долгота': 0.0
}

# Словарь переводов стран
country_translations = {
    'RU': 'Россия',
}

# Словарь переводов регионов
region_translations = {
    'Komi': 'Коми',
}

# Словарь переводов городов
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

# Сохраняем информацию о местоположении в переменную location_info
location_info = get_location_info(default_location)

# Выводим информацию о местоположении на экран
print(location_info)
