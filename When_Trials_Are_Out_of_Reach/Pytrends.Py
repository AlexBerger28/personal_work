from pytrends.request import TrendReq
import pandas as pd

# Define keywords
keywords = [
    "reiki", "energy healing", "ayahuasca", "sound bath",
    "acupuncture", "chiropractor", "naturopathy", "homeopathy",
    "herbal medicine", "functional medicine", "craniosacral therapy",
    "yoga retreat", "psychedelic integration", "holistic health",
    "meditation retreat"
]

# Connect to Google Trends
pytrends = TrendReq(hl='en-US', tz=360)

# Build payload
pytrends.build_payload(keywords, geo='US', timeframe='today 5-y')

# Get interest by metro area (DMA)
df = pytrends.interest_by_region(resolution='DMA', inc_low_vol=True, inc_geo_code=True)

# Sum total interest across keywords
df["total"] = df[keywords].sum(axis=1)

# Sort and save
df_sorted = df.sort_values(by="total", ascending=False)
df_sorted.to_csv("healing_hotspots.csv")
