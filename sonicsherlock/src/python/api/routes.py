"""
API routes for the audio fingerprinting service.
This file provides the web interface for the application.
"""

from flask import Blueprint, request, jsonify, current_app
import os
import tempfile
import librosa
import numpy as np
from werkzeug.utils import secure_filename
import time
import sys

# Add parent directory to path to import the main app
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from app import AudioFingerprinter

# Create blueprint
api = Blueprint('api', __name__)

# Initialize fingerprinter
fingerprinter = None

@api.before_app_first_request
def init_fingerprinter():
    """Initialize the fingerprinter on first request."""
    global fingerprinter
    db_path = current_app.config.get('DATABASE_PATH')
    fingerprinter = AudioFingerprinter(db_path=db_path)


@api.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        'status': 'ok',
        'message': 'Audio fingerprinting service is running'
    })


@api.route('/identify', methods=['POST'])
def identify():
    """Identify a song from uploaded audio file."""
    if 'file' not in request.files:
        return jsonify({
            'error': 'No file provided'
        }), 400
    
    file = request.files['file']
    
    if file.filename == '':
        return jsonify({
            'error': 'No file selected'
        }), 400
    
    # Save uploaded file to temporary location
    with tempfile.NamedTemporaryFile(delete=False, suffix='.wav') as temp:
        temp_path = temp.name
        file.save(temp_path)
    
    try:
        # Load audio
        start_time = time.time()
        y, sr = librosa.load(temp_path, sr=None, mono=True, duration=10.0)  # Limit to 10 seconds
        
        # Identify song
        result = fingerprinter.identify_song(y, sr)
        elapsed = time.time() - start_time
        
        if result:
            return jsonify({
                'matched': True,
                'title': result['title'],
                'artist': result['artist'],
                'album': result['album'],
                'confidence': round(result['confidence'], 2),
                'processing_time': round(elapsed, 2)
            })
        else:
            return jsonify({
                'matched': False,
                'message': 'No match found',
                'processing_time': round(elapsed, 2)
            })
    
    except Exception as e:
        return jsonify({
            'error': str(e)
        }), 500
    
    finally:
        # Clean up temporary file
        if os.path.exists(temp_path):
            os.unlink(temp_path)


@api.route('/add', methods=['POST'])
def add_song():
    """Add a song to the database."""
    if 'file' not in request.files:
        return jsonify({
            'error': 'No file provided'
        }), 400
    
    file = request.files['file']
    
    if file.filename == '':
        return jsonify({
            'error': 'No file selected'
        }), 400
    
    # Get metadata
    title = request.form.get('title')
    artist = request.form.get('artist')
    album = request.form.get('album')
    
    # If no title provided, use filename
    if not title:
        title = os.path.splitext(secure_filename(file.filename))[0]
    
    # Save uploaded file to temporary location
    with tempfile.NamedTemporaryFile(delete=False, suffix='.wav') as temp:
        temp_path = temp.name
        file.save(temp_path)
    
    try:
        # Add song to database
        start_time = time.time()
        song_id = fingerprinter.add_song(temp_path, title, artist, album)
        elapsed = time.time() - start_time
        
        return jsonify({
            'success': True,
            'song_id': song_id,
            'title': title,
            'artist': artist,
            'album': album,
            'processing_time': round(elapsed, 2)
        })
    
    except Exception as e:
        return jsonify({
            'error': str(e)
        }), 500
    
    finally:
        # Clean up temporary file
        if os.path.exists(temp_path):
            os.unlink(temp_path)


@api.route('/record', methods=['POST'])
def record():
    """Process audio recorded from the browser."""
    if 'audio_data' not in request.files:
        return jsonify({
            'error': 'No audio data provided'
        }), 400
    
    audio_file = request.files['audio_data']
    
    # Save uploaded file to temporary location
    with tempfile.NamedTemporaryFile(delete=False, suffix='.wav') as temp:
        temp_path = temp.name
        audio_file.save(temp_path)
    
    try:
        # Load audio
        y, sr = librosa.load(temp_path, sr=None, mono=True)
        
        # Identify song
        result = fingerprinter.identify_song(y, sr)
        
        if result:
            return jsonify({
                'matched': True,
                'title': result['title'],
                'artist': result['artist'],
                'album': result['album'],
                'confidence': round(result['confidence'], 2)
            })
        else:
            return jsonify({
                'matched': False,
                'message': 'No match found'
            })
    
    except Exception as e:
        return jsonify({
            'error': str(e)
        }), 500
    
    finally:
        # Clean up temporary file
        if os.path.exists(temp_path):
            os.unlink(temp_path)


@api.route('/stats', methods=['GET'])
def get_stats():
    """Get database statistics."""
    try:
        import sqlite3
        
        conn = sqlite3.connect(current_app.config.get('DATABASE_PATH'))
        cursor = conn.cursor()
        
        # Get song count
        cursor.execute("SELECT COUNT(*) FROM songs")
        song_count = cursor.fetchone()[0]
        
        # Get fingerprint count
        cursor.execute("SELECT COUNT(*) FROM fingerprints")
        fingerprint_count = cursor.fetchone()[0]
        
        # Get top 5 songs
        cursor.execute("""
            SELECT s.title, s.artist, COUNT(f.id) as fp_count 
            FROM songs s 
            JOIN fingerprints f ON s.id = f.song_id 
            GROUP BY s.id 
            ORDER BY fp_count DESC
            LIMIT 5
        """)
        top_songs = [{
            'title': row[0],
            'artist': row[1] or 'Unknown',
            'fingerprint_count': row[2]
        } for row in cursor.fetchall()]
        
        conn.close()
        
        return jsonify({
            'song_count': song_count,
            'fingerprint_count': fingerprint_count,
            'top_songs': top_songs
        })
    
    except Exception as e:
        return jsonify({
            'error': str(e)
        }), 500