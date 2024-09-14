# sentiment.py

from transformers import pipeline
import spacy

model1 = pipeline("sentiment-analysis", model="nlptown/bert-base-multilingual-uncased-sentiment")
nlp = spacy.load("ru_core_news_sm")

def get_sentiment(text, keywords):
    for sentiment, kw_list in keywords.items():
        if any(kw in text for kw in kw_list):
            return sentiment.capitalize()

    sentiment_result1 = model1(text)[0]
    label_map1 = {
        '1 star': 'Negative', '2 stars': 'Neutral', '3 stars': 'Neutral',
        '4 stars': 'Positive', '5 stars': 'Positive'
    }
    sentiment1 = label_map1.get(sentiment_result1['label'], 'Neutral')
    return sentiment1

def extract_location(text):
    doc = nlp(text)
    locations = [ent.text for ent in doc.ents if ent.label_ == "LOC"]
    return ', '.join(locations) if locations else None
