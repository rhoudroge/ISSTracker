utc = TimeScalesFactory.getUTC();
n = getCurrentTime(utc);
n = n.getComponents(utc).getDate;

fileName = [num2str(n.getYear) num2str(n.getMonth) num2str(n.getDay) '.txt']
