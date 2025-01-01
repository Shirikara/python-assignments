# This code generates a web app to calculate the area and circumference of circle and rectangle using Flask 

from flask import Flask, request, render_template_string
import math

app = Flask(__name__)

# HTML Template for user input
html_template = '''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Geometry Calculator</title>
    <style>
        body {
            background-color: #fcd6d6;
            font-family: Helvetica, sans-serif;
        }
        .image-container {
            position: absolute;
            top: 10px;
            right: 10px;
        }
        .image-container img {
            width: 200px;
            height: auto;
        }
    </style>
</head>
<body>
    <div class="image-container">
        <img src="https://andymath.com/area-and-circumference-circle-2/" alt="Circle Information">
    </div>
    <h1>Geometry Calculator</h1>

    <h2>Circle</h2>
    <form action="/calculate_circle" method="post">
        <label for="circle_radius">Radius:</label>
        <input type="number" step="0.01" name="circle_radius" id="circle_radius" required><br><br>
        <button type="submit">Calculate Circle</button>
    </form>

    <h2>Rectangle</h2>
    <form action="/calculate_rectangle" method="post">
        <label for="rec_width">Width:</label>
        <input type="number" step="0.01" name="rec_width" id="rec_width" required><br>
        <label for="rec_length">Length:</label>
        <input type="number" step="0.01" name="rec_length" id="rec_length" required><br><br>
        <button type="submit">Calculate Rectangle</button>
    </form>
    <!-- New image added here -->
    <div>
        <img src="https://andymath.com/wp-content/uploads/2019/08/circlenotes1.jpg" alt="Circle Notes">
         <img src="https://andymath.com/wp-content/uploads/2019/07/Rectangle-Notes-480x473.jpg" alt="Circle Notes">
    </div>
</body>
</html>
'''

# Route for input form
@app.route('/')
def index():
    return render_template_string(html_template)

# Route for circle calculations
@app.route('/calculate_circle', methods=['POST'])
def calculate_circle():
    try:
        circle_radius = float(request.form['circle_radius'])
        circle_area = math.pi * circle_radius ** 2
        circle_circumference = 2 * math.pi * circle_radius

        result = f'''
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <title>Circle Results</title>
            <style>
                body {{
                    background-color: lightpink;
                    font-family: Helvetica, sans-serif;
                }}
            </style>
        </head>
        <body>
            <h1>Circle Results</h1>
            <p>Radius: {circle_radius}</p>
            <p>Area: {round(circle_area, 2)}</p>
            <p>Circumference: {round(circle_circumference, 2)}</p>
            <br>
            <a href="/">Back to Calculator</a>
        </body>
        </html>
        '''

        return result

    except Exception as e:
        return f"<h1>Error</h1><p>{str(e)}</p><a href='/'>Back to Calculator</a>"

# Route for rectangle calculations
@app.route('/calculate_rectangle', methods=['POST'])
def calculate_rectangle():
    try:
        rec_width = float(request.form['rec_width'])
        rec_length = float(request.form['rec_length'])
        rec_area = rec_width * rec_length
        rec_circumference = 2 * (rec_width + rec_length)

        result = f'''
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <title>Rectangle Results</title>
            <style>
                body {{
                    background-color: lightpink;
                    font-family: Helvetica, sans-serif;
                }}
            </style>
        </head>
        <body>
            <h1>Rectangle Results</h1>
            <p>Width: {rec_width}</p>
            <p>Length: {rec_length}</p>
            <p>Area: {rec_area}</p>
            <p>Perimeter: {rec_circumference}</p>
            <br>
            <a href="/">Back to Calculator</a>
        </body>
        </html>
        '''

        return result

    except Exception as e:
        return f"<h1>Error</h1><p>{str(e)}</p><a href='/'>Back to Calculator</a>"

if __name__ == '__main__':
    app.run(debug=True)
