import geocoder

# Заглушки для местоположения
default_location = {
    'Страна': '',
    'Регион': '',
    'Город': '',
    'Широта': 0.0,
    'Долгота': 0.0
}

# Получаем информацию о местоположении
def get_location_info(defaults):
    try:
        g = geocoder.ip('me')

        if g.ok:
            return {
                'Страна': g.country if defaults['Страна'] == '' else defaults['Страна'],
                'Регион': g.state if defaults['Регион'] == '' else defaults['Регион'],
                'Город': g.city if defaults['Город'] == '' else defaults['Город'],
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
