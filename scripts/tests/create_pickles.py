import requests
from src.scrape import get_urls


def save_pickle(items):
    # Get first url for items
    urls = get_urls(items=items, n_pages=1)
    url = urls[0]

    # Save HTML and url for mock
    r = requests.get(url)
    out_file = f"tests/{items}"
    with open(f"{out_file}.html", "wb") as f:
        f.write(r.content)
    with open(f"{out_file}.txt", "w") as f:
        f.write(url)

    print(f"Saved: {out_file}.html")
    print(f"Saved: {out_file}.txt")


def main(items=None):
    print("Creating pickle files for testing.")
    if items is not None:
        save_pickle(items)
    else:
        items_list = ["author_urls", "paper_urls", "project_urls", "group_urls"]
        for items in items_list:
            save_pickle(items)


if __name__ == "__main__":
    main()
