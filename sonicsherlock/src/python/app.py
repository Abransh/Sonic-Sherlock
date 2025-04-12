#!/usr/bin/env python3
"""
Main application for audio fingerprinting service.
This combines Python ease-of-use with C++ performance.
"""

import os
import sys
import numpy as np
import time
import argparse
import sqlite3
from pathlib import Path

# Import C++ extension module (will be compiled with pybind11)
try:
    import audiofp_core as afp
except ImportError:
    print("C++ extension module not found. Please build it first with:")
    print("    python setup.py build_ext --inplace")
    sys.exit(1)

# For audio file handling
import librosa


class AudioFingerprinter:
    """Main class for audio fingerprinting functionality."""
    
    def __init__(self, db_path=None):
        """Initialize the fingerprinter with optional database path."""
        self.db_path = db_path or os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'db', 'fingerprints.db')
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        self._init_db()
    
    def _init_db(self):
        """Initialize the database if it doesn't exist."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Create tables if they don't exist
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS songs (
            id INTEGER PRIMARY KEY,
            title TEXT NOT NULL,
            artist TEXT,
            album TEXT,
            duration REAL
        )
        ''')
        
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS fingerprints (
            id INTEGER PRIMARY KEY,
            song_id INTEGER,
            hash INTEGER NOT NULL,
            time_offset INTEGER NOT NULL,
            FOREIGN KEY (song_id) REFERENCES songs(id)
        )
        ''')
        
        # Create index for faster lookups
        cursor.execute('''
        CREATE INDEX IF NOT EXISTS idx_fingerprints_hash ON fingerprints(hash)
        ''')
        
        conn.commit()
        conn.close()
    
    def fingerprint_file(self, file_path):
        """Generate fingerprints from an audio file using C++ for performance."""
        print(f"Fingerprinting {file_path}...")
        start_time = time.time()
        
        # Load audio file using librosa (keep in Python for convenience)
        y, sr = librosa.load(file_path, sr=None, mono=True)
        
        # Convert to numpy array of correct type for C++
        samples = np.array(y, dtype=np.float64)
        
        # Use C++ for performance-critical operations
        print("Generating spectrogram...")
        spectrogram = afp.Fingerprinter.generate_spectrogram(samples, sr)
        
        print("Finding peaks...")
        peaks = afp.PeakFinder.find_peaks(spectrogram, neighborhood_size=10, threshold=0.3, max_peaks=100)
        
        print("Generating fingerprints...")
        fingerprints = afp.Fingerprinter.generate_fingerprints(spectrogram, sr)
        
        elapsed = time.time() - start_time
        print(f"Generated {len(fingerprints)} fingerprints in {elapsed:.2f} seconds")
        
        return fingerprints
    
    def add_song(self, file_path, title=None, artist=None, album=None):
        """Add a song to the database."""
        # Extract metadata if not provided
        if title is None:
            title = Path(file_path).stem
        
        # Get audio duration
        y, sr = librosa.load(file_path, sr=None, duration=30)  # Just load a bit to get duration
        duration = librosa.get_duration(y=y, sr=sr)
        
        # Insert song into database
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute(
            "INSERT INTO songs (title, artist, album, duration) VALUES (?, ?, ?, ?)",
            (title, artist, album, duration)
        )
        song_id = cursor.lastrowid
        
        # Generate fingerprints
        fingerprints = self.fingerprint_file(file_path)
        
        # Batch insert fingerprints for performance
        print(f"Inserting {len(fingerprints)} fingerprints into database...")
        batch_size = 1000
        for i in range(0, len(fingerprints), batch_size):
            batch = fingerprints[i:i+batch_size]
            values = [(song_id, fp.hash, fp.time_offset) for fp in batch]
            cursor.executemany(
                "INSERT INTO fingerprints (song_id, hash, time_offset) VALUES (?, ?, ?)",
                values
            )
        
        conn.commit()
        conn.close()
        
        print(f"Added song '{title}' with ID {song_id} to database")
        return song_id
    
    def identify_song(self, audio_samples, sample_rate):
        """Identify a song from audio samples."""
        print("Generating fingerprints from audio sample...")
        
        # Convert to numpy array
        samples = np.array(audio_samples, dtype=np.float64)
        
        # Generate fingerprints using C++ code
        fingerprints = afp.Fingerprinter.fingerprint_audio(samples, sample_rate)
        
        print(f"Generated {len(fingerprints)} fingerprints from sample")
        
        # Search for matches in database
        matches = self._find_matches(fingerprints)
        
        if not matches:
            print("No matches found")
            return None
        
        # Get the top match
        top_match = matches[0]
        print(f"Top match: {top_match['title']} by {top_match['artist']} ({top_match['confidence']:.2f}% confidence)")
        
        return top_match
    
    def _find_matches(self, fingerprints):
        """Find matches for a set of fingerprints in the database."""
        if not fingerprints:
            return []
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Find all matching fingerprints
        matches = {}
        for fp in fingerprints:
            cursor.execute(
                "SELECT song_id, time_offset FROM fingerprints WHERE hash = ?",
                (fp.hash,)
            )
            
            for song_id, db_offset in cursor.fetchall():
                time_delta = db_offset - fp.time_offset
                
                if song_id not in matches:
                    matches[song_id] = {}
                
                if time_delta not in matches[song_id]:
                    matches[song_id][time_delta] = 0
                
                matches[song_id][time_delta] += 1
        
        # Calculate best matching offset for each song
        results = []
        for song_id, offsets in matches.items():
            # Find offset with most matches
            max_matches = max(offsets.values())
            
            # Get song details
            cursor.execute(
                "SELECT title, artist, album FROM songs WHERE id = ?",
                (song_id,)
            )
            song_info = cursor.fetchone()
            
            # Calculate confidence score (simple version)
            confidence = (max_matches / len(fingerprints)) * 100
            
            # Add to results
            results.append({
                'song_id': song_id,
                'title': song_info[0],
                'artist': song_info[1] or 'Unknown',
                'album': song_info[2] or 'Unknown',
                'matches': max_matches,
                'confidence': confidence
            })
        
        conn.close()
        
        # Sort by number of matches
        results.sort(key=lambda x: x['matches'], reverse=True)
        
        return results


def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(description='Audio Fingerprinting Tool')
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Add command
    add_parser = subparsers.add_parser('add', help='Add a song to the database')
    add_parser.add_argument('file_path', help='Path to audio file')
    add_parser.add_argument('--title', help='Song title')
    add_parser.add_argument('--artist', help='Artist name')
    add_parser.add_argument('--album', help='Album name')
    
    # Identify command
    identify_parser = subparsers.add_parser('identify', help='Identify a song from audio file')
    identify_parser.add_argument('file_path', help='Path to audio file to identify')
    identify_parser.add_argument('--duration', type=float, default=10.0, help='Duration in seconds to sample')
    
    args = parser.parse_args()
    
    fingerprinter = AudioFingerprinter()
    
    if args.command == 'add':
        fingerprinter.add_song(args.file_path, args.title, args.artist, args.album)
    
    elif args.command == 'identify':
        # Load audio file
        y, sr = librosa.load(args.file_path, sr=None, mono=True, duration=args.duration)
        
        # Identify
        result = fingerprinter.identify_song(y, sr)
        
        if result:
            print(f"Found: {result['title']} by {result['artist']}")
            print(f"Confidence: {result['confidence']:.2f}%")
        else:
            print("No match found")
    
    else:
        parser.print_help()


if __name__ == '__main__':
    main()