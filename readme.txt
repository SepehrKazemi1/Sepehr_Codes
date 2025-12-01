# Drought–Heatwave Analysis (Central Iran)

This repo contains my MSc work on compound drought–heatwave events in the Zayandehrud Basin and Isfahan province (1980–2024).

## What this code does

- Reads Era5 data and station time series
- Calculates SPI and SPEI at different time scales
- Derives several heatwave indices (e.g. duration, frequency, intensity)
- Makes basic plots and maps for checking results

The scripts are written in MATLAB. Some helper tools in Python are used for cleaning data and exporting tables.

## Why I wrote this

I needed a workflow that I can re run when I fix data issues or change thresholds. A lot of the code is “research style”: it works, but it’s not a polished package.

## How to use

1. Put your climate data in `data/` (see `data/README.md` for format).
2. Edit `config.m` with your paths and basin settings.
3. Run `main_spi_spei.m` and then `main_heatwaves.m`.

## Limitations

- Hard-coded for my study area (central Iran).
- Error handling is minimal.
- Some variable names are in Persian because that’s how I first wrote them 🙂

I’ll clean this up more as I go through my thesis.
