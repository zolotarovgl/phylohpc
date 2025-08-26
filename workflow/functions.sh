# Define the function to read and echo configuration variables
read_config() {
    local config_file="$1"

    if [[ ! -f "$config_file" ]]; then
        echo "Error: Config file '$config_file' not found!"
        return 1
    fi

    #echo "------------------"
    #echo -e "CONFIG_FILE: ${config_file}\n"

    # Read variables from the configuration file
    while IFS='=' read -r key value; do
        # Skip empty lines and comments
        if [[ -n "$key" && ! "$key" =~ ^# ]]; then
            eval "${key}='${value}'"
            # Echo the variable and its value
            #echo "$key=$value"
        fi
    done < "$config_file"

    #echo "------------------"
}


check_file_exists() {
    local file="$1"

    if [[ ! -f "$file" ]]; then
        echo "Error: File '$file' does not exist."
        exit 1
    else
        echo "$file exists"
    fi
}

get_ram() {
    sacct --format=JobID,MaxRSS%9 -j "$1"
}

