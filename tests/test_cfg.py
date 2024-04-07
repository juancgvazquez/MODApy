import os
import unittest
from unittest.mock import patch

from MODApy.cfg import Config


class TestConfig(unittest.TestCase):
    @patch('MODApy.cfg.Redis')
    def test_redis_connection(self, mock_redis):
        # Configure the mock Redis
        mock_redis_instance = mock_redis.return_value
        config = Config()
        # Check if Redis connection is correctly initialized
        self.assertEqual(config.redis_conn, mock_redis_instance)

    def test_set_config(self):
        config = Config()
        # Test setting a configuration value
        config.setConfig("GENERAL", "testmode", "False")
        self.assertEqual(config.cfg["GENERAL"]["testmode"], "False")

    def test_default_values(self):
        config = Config()
        expected_root_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'MODApy'
        )
        # Test default values
        self.assertEqual(config.rootDir, expected_root_dir)
        self.assertEqual(config.testFlag, False)

    @patch('MODApy.cfg.logging')
    def test_logging_setup(self, mock_logging):
        config = Config()
        config.setup_logging()
        # Check if logging configurations are set correctly
        mock_logging.config.dictConfig.assert_called_once()


if __name__ == '__main__':
    unittest.main()
