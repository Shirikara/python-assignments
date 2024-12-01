# ------------------- This is a guessing a number user interactive game ------------------- #

# Import Libraries
import numpy as np

# Define the game inside a function
def guess_number_game():
    print("Welcome to the telepathy game!")
    # Initializing Game Variables
    total_num_games = 0 #keeps track of how many games the user has played
    play_again = True #controls the main game. It starts as True so the user can play at least one game

    while play_again:
        # Generate a random number for the game
        comp_int_1to20 = np.random.randint(1, 21) #random int between 1 and 20 
        total_num_games += 1 #counts if new game is starting
        num_try = 0 #How many times it takes the user to guess correctly within one round
        keep_playing = True #This loop keeps running as long as the user wants to play (keep_playing = True)

        #Game instructions and special commands for the user
        print("\nI am thinking of a number between 1 and 20. Can you guess what it is?")
        print("Enter 'x' to exit the program, 'n' to start a new game, or 's' to show the hidden number.")

        while keep_playing: #This loop continues until the current game ends (keep_playing = False).
            user_input = input("Your guess: ")
            # Special Commands
            if user_input.lower() == 'x':
                print("Thanks for playing with me. Goodbye!")
                return
            elif user_input.lower() == 'n':
                print("Starting a new game!")
                keep_playing = False
            elif user_input.lower() == 's':
                print(f"Hey, you are cheating! The hidden number is {comp_int_1to20}")
            else:
                try:
                    guess_by_user = int(user_input) # If the input is not a command (x, n, s), the program attempts to convert it to an integer
                    num_try += 1
                    
                    #Feedback to the user's guess
                    if guess_by_user == comp_int_1to20:
                        print(f"Yay! You guessed correctly in {num_try} tries. You read my thoughts and won the game!")
                        keep_playing = False
                    elif guess_by_user < 0:
                        print("Only positive values are accepted. Please try again")
                    elif guess_by_user == 0:
                        print("Invalid Input. The minimal value is 1. Please try again")
                    elif guess_by_user > 20:
                        print("Invalid Input. The maximal value is 20. Please try again")
                    elif guess_by_user < comp_int_1to20:
                        print("Your number is too small. Try again.")
                    else:
                        print("Your number is too big. Try again.")
                except ValueError:
                    print("Invalid input. Please enter a positive integer (natural number, no fractions) between 1 and 20, or 'x', 'n', 's'.")

        # Ask if the user wants to play again after the current game ends
        if keep_playing == False:
            play_again_input = input("Would you like to play again? (yes/no): ").strip().lower()
            if play_again_input not in ['yes', 'y']:
                play_again = False
    #Exit the Game
    print(f"Thanks for playing! You played {total_num_games} games. Have a nice day :)")

#Initialize the game by running the function
guess_number_game()