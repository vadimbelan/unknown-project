# scraper.py

import requests
from bs4 import BeautifulSoup
from .db import add_news_to_db
from .sentiment import get_sentiment, extract_location

def parse_news(source, latitude, longitude, table, search_url, news_selector, title_selector, keywords):
    response = requests.get(search_url.format(source=source))
    soup = BeautifulSoup(response.content, 'html.parser')

    news_list = []
    for item in soup.select(news_selector)[:10]:
        title = item.select_one(title_selector).get_text(strip=True) if title_selector else item.get_text(strip=True)
        sentiment = get_sentiment(title, keywords)
        location = extract_location(title)
        add_news_to_db(table, title, latitude, longitude, sentiment, location)
        news_list.append({
            'title': title, 'latitude': latitude, 'longitude': longitude, 'sentiment': sentiment, 'location': location
        })
    return news_list
