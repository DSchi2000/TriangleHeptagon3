import smtplib
import numpy as np


def sendMail(message):

    YOUR_GOOGLE_EMAIL = 'themaker845@gmail.com'  
    YOUR_GOOGLE_EMAIL_APP_PASSWORD = 'klja zprt fodh phcy'  

    smtpserver = smtplib.SMTP_SSL('smtp.gmail.com', 465)
    smtpserver.ehlo()
    smtpserver.login(YOUR_GOOGLE_EMAIL, YOUR_GOOGLE_EMAIL_APP_PASSWORD)

    sent_from = YOUR_GOOGLE_EMAIL
    sent_to = "daniel.schieber@ugent.be"  
    SUBJECT='[LANAYRU] CODE FINISHED'
    TEXT=message
    email_text = f'Subject: {SUBJECT}\n\n{TEXT}'
    smtpserver.sendmail(sent_from, sent_to, email_text)

    smtpserver.close()

timeused=np.loadtxt("TIME.txt",dtype="str",delimiter=";")
message=f"RUN FINISHED\n----{timeused[0]}---- \nTIME USED: {float(timeused[1]):.0f} s"
sendMail(message)