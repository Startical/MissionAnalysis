from datetime import datetime, timedelta


def calendar_date_to_datetime(date_string):
    for fmt in ('%Y-%m-%dT%H:%M:%S.%fZ', '%Y-%m-%dT%H:%M:%SZ'):
        try:
            return datetime.strptime(date_string, fmt)
        except ValueError:
            continue
        raise ValueError(f"Time data '{date_string}' does not match expected formats.")

def time_offset(date1_string, date2_string):

    date1 = calendar_date_to_datetime(date1_string)
    date2 = calendar_date_to_datetime(date2_string)

    return (date1-date2).total_seconds()



