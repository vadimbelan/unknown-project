import geocoder

def get_location_info():
    try:
        g = geocoder.ip('me')
        if g.ok:
            return {
                'IP': g.ip,
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

location_info = get_location_info()
print(location_info)
