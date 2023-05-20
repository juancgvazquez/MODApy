import logging
import os
import unittest

from MODApy.utils import InvalidFileError, checkFile

logger = logging.getLogger(__name__)


class TestCheckFile(unittest.TestCase):
    def test_existing_file_with_correct_extension(self):
        # Arrange
        file_path = "test_file.txt"
        file_extension = ".txt"
        with open(file_path, "w") as file:
            file.write("Test content")

        # Act
        result = checkFile(file_path, file_extension)

        # Assert
        self.assertTrue(result)

        # Cleanup
        os.remove(file_path)

    def test_existing_file_with_incorrect_extension(self):
        # Arrange
        file_path = "test_file.txt"
        file_extension = ".pdf"
        with open(file_path, "w") as file:
            file.write("Test content")

        # Act and Assert
        with self.assertRaises(InvalidFileError):
            checkFile(file_path, file_extension)

        # Cleanup
        os.remove(file_path)

    def test_non_existing_file(self):
        # Arrange
        file_path = "non_existing_file.txt"
        file_extension = ".txt"

        # Act and Assert
        with self.assertRaises(InvalidFileError):
            checkFile(file_path, file_extension)


if __name__ == '__main__':
    unittest.main()
