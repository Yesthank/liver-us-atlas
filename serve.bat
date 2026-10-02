@echo off
rem Start a local web server for the atlas and open it in the default browser.
cd /d "%~dp0"
start "" http://localhost:8917/
python -m http.server 8917
