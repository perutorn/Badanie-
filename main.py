from application import appGUI
import multiprocessing
import settings
import os
import pathlib


if __name__ == "__main__":
    #freeze_support jest niezbędne aby program mógł zostać spakowany jako
    #EXE pod Windowsem. Nie ma znaczenia na innych platformach
    multiprocessing.freeze_support()

    #wygeneruj plik ustawień domyślnych eżeli nie znaleziono żadnegi innego pliku ustawień

    settings_file = 'settings.json'
    path = pathlib.Path(__file__).parent.absolute()

    path = str(path / settings_file)
    if os.path.exists(path):
        print("Znaleziono settings.json")
    else:
        print("Nieznaleziono settings.json: Generowanie pliku domyślengo")
        settings.Settings.from_defaults()
    #punkt wejścia do całęgo programu:
    app = appGUI()
    app.run()