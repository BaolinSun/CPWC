import logging


def create_logger(infoname):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s %(levelname)s PROCESS#%(process)d %(message)s')

    console_handler = logging.StreamHandler()
    file_handler = logging.FileHandler(filename=infoname, encoding='UTF-8')
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger