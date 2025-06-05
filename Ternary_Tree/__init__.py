import logging.config
from .ucc.upgccsd import CircWrapper, LadExcImpl
import logging

logging.FileHandler()
DEFAULT_LOGGING_CONFIG = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'standard': {'format': '%(asctime)s [%(levelname)s]: %(message)s'}
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'standard',
            'level': 'INFO'
        },
        'file': {
            'class': 'logging.FileHandler',
            'filename': f'{__name__}/package.log',
            'formatter': 'standard',
            'level': 'DEBUG',
        }
    },
    'loggers': {
        f'{__name__}': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': True
        },
    }
}

logging.config.dictConfig(DEFAULT_LOGGING_CONFIG)

logger = logging.getLogger(__name__) 