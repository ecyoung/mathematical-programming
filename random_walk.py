"""
random_walk.py
"""

import math # import the math module
import random # import the random module

# problem 1
def rps(): # define function with no inputs that returns nothing

	rpsList = ["rock", "paper", "scissors"] # create a list of possible outputs
	gameCounter = 0 # initialize game counter
	winCounter = 0 # initialize win counter

	print('A challenger appears!')
	intro = input("Enter your name: ") # prompts user to input name
	print("Welcome " + intro + "!") # welcomes user

	while True: 
		# y = yes, n = no
		contStatus = input("Enter y or n to indicate whether or not you want to play a game of rock-paper-scissors: ")
		# if the user selects yes, increment game counter by 1, prompt user to enter output, randomly select output for opponent, and print opponenet output
		if contStatus == "y":
			gameCounter += 1
			userChoice = input("Enter rock, paper, or scissors: ")
			oppChoice = random.choice(rpsList)
			print("Your opponent used", oppChoice)
			# for all possible combinations, print the game result, number of games played, and number of games won
			# note that rock beats scissors, scissors beats paper, and paper beats rock
			if userChoice == "rock" and oppChoice == "rock":
				print("It's a tie!")
				print("Games Played:", gameCounter)
				print("Games Won:", winCounter)

			elif userChoice == "rock" and oppChoice == "paper":
				print("You lose :(")
				print("Games Played:", gameCounter)
				print("Games Won:", winCounter)

			elif userChoice == "rock" and oppChoice == "scissors":
				print("You win :)")
				winCounter += 1 # increment win counter by 1 whenever user wins
				print("Games Played:", gameCounter)
				print("Games Won:", winCounter)

			elif userChoice == "paper" and  oppChoice == "rock":
				print("You win :)")
				winCounter += 1
				print("Games Played:", gameCounter)
				print("Games Won:", winCounter)

			elif userChoice == "paper" and oppChoice == "paper":
				print("It's a tie!")
				print("Games Played:", gameCounter)
				print("Games Won:", winCounter)

			elif userChoice == "paper" and oppChoice == "scissors":
				print("You lose :(")
				print("Games Played:", gameCounter)
				print("Games Won:", winCounter)

			elif userChoice == "scissors" and oppChoice == "rock":
				print("You lose :(")
				print("Games Played:", gameCounter)
				print("Games Won:", winCounter)

			elif userChoice == "scissors" and oppChoice == "paper":
				print("You win!")
				winCounter += 1
				print("Games Played:", gameCounter)
				print("Games Won:", winCounter)

			elif userChoice == "scissors" and oppChoice == "scissors":
				print("It's a tie!")
				print("Games Played:", gameCounter)
				print("Games Won:", winCounter)
		# if user selects to not play another game, print the game history and break from the while loop
		if contStatus == "n": 
			print('Goodbye, it was fun playing!')
			print("Games Played:", gameCounter)
			print("Games Won:", winCounter)
			break 
	return

# problem 2
def sleepyStudent(leftProb, dist): # leftProb = probability of moving left 1 unit, dist = distance to either end from center
    numSteps = 0
    location = dist
 
    while location != 0 or location != 2*dist: # loop while the location is not at the start or end of the walkway
        threshold = random.random()
        if threshold <= leftProb: # random() selects a random float between 0 to 1
            location -= 1 # move 1 unit to the left
            numSteps += 1 # increment number of steps taken by 1
        else:
            location += 1 # move 1 unit to the right
            numSteps += 1 
        walkway = '-'*location + 'S' + '-'*(2*dist-location) # S represents the student
        print('|' + walkway + '|') # illustrate current position in walkway
        
        if location == 0 or location == 2*dist: # if the student is at the start or end of the walkway
            print('Student took %d steps.' % numSteps) # print the number of steps taken
            break # break from while loop

    return numSteps # returns number of steps taken