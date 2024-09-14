import psycopg2

def create_user_table(cur):
    cur.execute("DROP TABLE IF EXISTS users")
    cur.execute('''CREATE TABLE IF NOT EXISTS users (
                    id SERIAL PRIMARY KEY,
                    username TEXT UNIQUE NOT NULL,
                    email TEXT UNIQUE NOT NULL,
                    password TEXT NOT NULL,
                    budget REAL DEFAULT NULL,
                    currency TEXT DEFAULT 'руб.'
                  )''')

def add_user(username, email, password):
    conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")
    with conn:
        cur = conn.cursor()
        cur.execute("INSERT INTO users (username, email, password) VALUES (%s, %s, %s)",
                    (username, email, password))
    conn.commit()

def get_user_by_email(email):
    conn = psycopg2.connect("dbname=adventNews user=postgres password=1234 host=localhost")
    with conn:
        cur = conn.cursor()
        cur.execute("SELECT * FROM users WHERE email = %s", (email,))
        return cur.fetchone()
