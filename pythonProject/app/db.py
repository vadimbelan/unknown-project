#db.py

import psycopg2
from .models import create_user_table

def init_db():
    conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")
    with conn:
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS news_nova_rambler_ru CASCADE")
        cur.execute('''CREATE TABLE news_nova_rambler_ru (
                        id SERIAL PRIMARY KEY,
                        title TEXT,
                        latitude REAL,
                        longitude REAL,
                        sentiment TEXT,
                        location TEXT)''')
        cur.execute("DROP TABLE IF EXISTS news_tags CASCADE")
        cur.execute('''CREATE TABLE news_tags (
                        id SERIAL PRIMARY KEY,
                        title TEXT,
                        latitude REAL,
                        longitude REAL,
                        sentiment TEXT,
                        location TEXT)''')

        create_user_table(cur)
        create_favorites_table(cur)

    conn.commit()

def add_news_to_db(table, title, latitude, longitude, sentiment, location):
    conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")
    with conn:
        cur = conn.cursor()
        cur.execute(f"INSERT INTO {table} (title, latitude, longitude, sentiment, location) VALUES (%s, %s, %s, %s, %s)",
                    (title, latitude, longitude, sentiment, location))
        cur.execute(f"DELETE FROM {table} WHERE id NOT IN (SELECT id FROM {table} ORDER BY id DESC LIMIT 10)")
    conn.commit()

def create_favorites_table(cur):
    cur.execute("DROP TABLE IF EXISTS favorites CASCADE")
    cur.execute('''CREATE TABLE IF NOT EXISTS favorites (
                        id SERIAL PRIMARY KEY,
                        news_id INTEGER,
                        FOREIGN KEY(news_id) REFERENCES news_nova_rambler_ru(id))''')
