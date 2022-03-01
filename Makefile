batch := ""

.PHONY: script
script:
	@[ "${nb}" ] && jupyter nbconvert --to script $(nb) || ( echo "nb variable not set."; return 1 )

scrape:
	@[ "${items}" ] && python 0_async_scrape.py $(items) $(batch)|| (echo "target variable not set."; return 1 )
