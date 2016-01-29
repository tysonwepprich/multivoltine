# email me when analysis finishes!!!

list.of.packages <- "mailR"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)
library(mailR)

sess <- sessionInfo()
send.mail(from = "tysonR@gmail.com",
          to = "Recipient 1 <tyson.wepprich@gmail.com>",
          subject = "This is your R speaking",
          body = paste("Your analysis on", sess$running, "is finished", sep = " "),
          smtp = list(host.name = "aspmx.l.google.com", port = 25),
          authenticate = FALSE,
          send = TRUE)