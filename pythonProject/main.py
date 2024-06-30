import geocoder


# Получаем информацию о местоположении
def get_location_info():
    try:
        g = geocoder.ip('me')

        if g.ok:
            return {
                'Страна': g.country,
                'Регион': g.state,
                'Город': g.city,
                'Широта': g.latlng[0],
                'Долгота': g.latlng[1]
            }
        else:
            return {'Ошибка': 'Не удалось получить местоположение'}
    except Exception as e:
        return {'Ошибка': str(e)}


# Сохраняем информацию о местоположени в переменную location_info
location_info = get_location_info()

# Выводим информацию о местоположении на экран
print(location_info)
