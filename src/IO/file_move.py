import os 
import shutil

def andes_report_transfer(dir2):
    target_folder = os.path.join(dir2, "andes_report")

    # Specify the target folder where .txt files should be moved
    # target_folder = r"C:\Users\xxx\Documents\desired_folder"  # Change to your desired path

    # Create the target folder if it doesn't exist
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)
        print(f"Folder created: {target_folder}")
    else:
        print(f"Folder already exists: {target_folder}")

    # Get the current working directory
    current_dir = os.getcwd() #os.path.dirname(dir1)

    # Iterate over all files in the current working directory
    for file_name in os.listdir(current_dir):
        # Check if the file is a .txt file
        if file_name.endswith(".txt"):
            source_path = os.path.join(current_dir, file_name)  # Full path to the file
            target_path = os.path.join(target_folder, file_name)  # Full path in the target folder
            # Move the .txt file
            shutil.move(source_path, target_path)
            print(f"Moved file: {file_name} -> {target_path}")

        if file_name.endswith(".npz"):
            source_path = os.path.join(current_dir, file_name)  # Full path to the file
            target_path = os.path.join(target_folder, file_name)  # Full path in the target folder
            # Move the .npz file
            shutil.move(source_path, target_path)
            print(f"Moved file: {file_name} -> {target_path}")

        if file_name.endswith(".lst"):
            source_path = os.path.join(current_dir, file_name)  # Full path to the file
            target_path = os.path.join(target_folder, file_name)  # Full path in the target folder

            # Move the .lst file
            shutil.move(source_path, target_path)
            print(f"Moved file: {file_name} -> {target_path}")