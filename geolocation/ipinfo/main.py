import requests

def get_location_info():
    try:
        response = requests.get('https://ipinfo.io/json')
        data = response.json()

        # Проверяем наличие необходимых полей в ответе
        if 'ip' in data:
            ip = data.get('ip', 'Неизвестно')
            city = data.get('city', 'Неизвестно')
            region = data.get('region', 'Неизвестно')
            country = data.get('country', 'Неизвестно')
            loc = data.get('loc', '0,0').split(',')
            latitude = loc[0]
            longitude = loc[1]

            return {
                'IP': ip,
                'Страна': country,
                'Регион': region,
                'Город': city,
                'Широта': latitude,
                'Долгота': longitude
            }
        else:
            return {'Ошибка': 'Не удалось получить местоположение'}
    except Exception as e:
        return {'Ошибка': str(e)}

location_info = get_location_info()
print(location_info)
