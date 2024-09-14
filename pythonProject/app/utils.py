#utils.py

import requests
import psycopg2
import folium
from flask import render_template

CIRCLE_RADIUS = 50000

def get_location_info(defaults):
    try:
        response = requests.get("https://ipinfo.io/json")
        data = response.json()
        loc = data.get('loc', '0,0').split(',')
        return {
            'Страна': defaults['Страна'] if defaults['Страна'] else data.get('country', ''),
            'Регион': defaults['Регион'] if defaults['Регион'] else data.get('region', ''),
            'Город': defaults['Город'] if defaults['Город'] else data.get('city', ''),
            'Широта': defaults['Широта'] if defaults['Широта'] else float(loc[0]),
            'Долгота': defaults['Долгота'] if defaults['Долгота'] else float(loc[1])
        }
    except Exception as e:
        defaults['Ошибка'] = str(e)
    return defaults



def write_html_for_two_cities(content, loc_info_user, loc_info_geo, src, ts, tags,
                              circle_color_user, sentiment_message_user,
                              circle_color_geo, sentiment_message_geo):
    # Создаем карту с центром на среднем значении между двумя городами
    map_center = [(loc_info_user['Широта'] + loc_info_geo['Широта']) / 2, (loc_info_user['Долгота'] + loc_info_geo['Долгота']) / 2]
    my_map = folium.Map(location=map_center, zoom_start=5)

    # Добавляем круг для введенного города
    folium.Circle(
        location=[loc_info_user['Широта'], loc_info_user['Долгота']],
        radius=CIRCLE_RADIUS,
        color=circle_color_user,
        fill=True,
        fill_color=circle_color_user,
        fill_opacity=0.3,
        stroke=True,
        opacity=0.5
    ).add_to(my_map)

    # Добавляем круг для города по геопозиции
    folium.Circle(
        location=[loc_info_geo['Широта'], loc_info_geo['Долгота']],
        radius=CIRCLE_RADIUS,
        color=circle_color_geo,
        fill=True,
        fill_color=circle_color_geo,
        fill_opacity=0.3,
        stroke=True,
        opacity=0.5
    ).add_to(my_map)

    # Рассчитываем среднюю точку между двумя городами
    midpoint = [
        (loc_info_user['Широта'] + loc_info_geo['Широта']) / 2,
        (loc_info_user['Долгота'] + loc_info_geo['Долгота']) / 2
    ]

    # Добавляем линии на карте
    folium.PolyLine(
        locations=[[loc_info_user['Широта'], loc_info_user['Долгота']], midpoint],
        color=circle_color_user,
        weight=2.5,
        opacity=1
    ).add_to(my_map)

    folium.PolyLine(
        locations=[midpoint, [loc_info_geo['Широта'], loc_info_geo['Долгота']]],
        color=circle_color_geo,
        weight=2.5,
        opacity=1
    ).add_to(my_map)

    # Сохраняем карту в HTML
    map_html = my_map._repr_html_()

    return render_template('two_cities.html', src=src, ts=ts, tags=tags, content=content,
                           sentiment_message_user=sentiment_message_user,
                           sentiment_message_geo=sentiment_message_geo,
                           map_html=map_html)


from flask import render_template

def view_news_table(table_name, title):

    conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")

    with conn:
        cur = conn.cursor()
        cur.execute(f"SELECT * FROM {table_name}")
        rows = cur.fetchall()

    return render_template('news_table.html', title=title, rows=rows)


def is_news_favorite(news_id):

    conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")

    with conn:
        cur = conn.cursor()
        cur.execute("SELECT 1 FROM favorites WHERE news_id = %s", (news_id,))
        return cur.fetchone() is not None

keywords = {
    'positive': ['Пожар потушили жители вместе со спасателями в Автозаводском районе'],
    'neutral': ['Маслкар-легенда Chevrolet Camaro сбил девушку на зебре в Новосибирске'],
    'negative': ['Синоптик: Сентябрь в Москве может оказаться рекордно сухим']
}
