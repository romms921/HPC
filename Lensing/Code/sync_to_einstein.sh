#!/bin/bash

# --- Configuration ---
# Set the path where the IP address file will be saved.
IP_FILE_PATH="/Users/ainsleylewis/Documents/Astronomy/HPC/Lensing/Code/public_ip.txt"

# Set the local source directory you want to sync.
LOCAL_SOURCE_DIR="/Users/ainsleylewis/Documents/Astronomy/HPC/Lensing/"

# Set the remote destination on your server.
REMOTE_DESTINATION="rommulus@hpc2021-io1.hku.hk:/home/rommulus/Projects/itng_lensing"

# --- Main Script Logic ---

echo "🔎 Step 1: Finding your public IP address..."

# Use an online service to get the public IP. The '-s' flag makes curl silent.
PUBLIC_IP=$(ifconfig utun6 | grep 'inet ' | awk '{print $2}')

# Check if curl successfully retrieved an IP address.
if [ -z "$PUBLIC_IP" ]; then
    echo "❌ Error: Could not retrieve public IP address. Please check your internet connection."
    exit 1
fi

echo "✅ Your public IP is: $PUBLIC_IP"

# Ensure the directory for the IP file exists.
mkdir -p "$(dirname "$IP_FILE_PATH")"

echo "💾 Step 2: Saving IP to $IP_FILE_PATH..."
echo "$PUBLIC_IP" > "$IP_FILE_PATH"

echo "🚀 Step 3: Starting rsync to send files to Einstein..."
rsync -avzP "$LOCAL_SOURCE_DIR" "$REMOTE_DESTINATION"

echo "🎉 Sync to Einstein complete."