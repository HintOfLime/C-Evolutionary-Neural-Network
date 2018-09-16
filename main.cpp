#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <SFML/Graphics.hpp>

const int WIDTH = 800;
const int HEIGHT = 400;

const int COLUMNS = 4;
const int ROWS = 15;

const int POPULATION = 50;
const int TRIALS = 50;

const int OFFSPRING_PER_GENERATION = 30;
const int MUTATIONS = 5;

struct Neuron {
        Neuron* output[ROWS];
        float weight[ROWS];
        float value;
};

Neuron* neurons [COLUMNS][ROWS];
float genes[POPULATION][COLUMNS][ROWS][ROWS];

float targetFunction (float input) {
    //float theta = (input+1)*M_PI;
    //return sin(theta);
    return input;
}

float activationFunction (float input) {
    return std::max(tanh(input), (double)0);
}

void swap(int *xp, int *yp) {
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void swap(float *xp, float *yp) {
    float temp = *xp;
    *xp = *yp;
    *yp = temp;
}
 

float simulate (float input) {
    // Reset neuron values
    for (int i = 0; i < COLUMNS; i++) {
        for (int j = 0; j < ROWS; j++) {
            for (int k = 0; k < ROWS; k++) {
                neurons[i][j]->value = 0;
            }
        }
    }

    // Set input layer values
    neurons[0][0]->value = input;
    neurons[0][ROWS-1]->value = 1;

    // Simulate network
    for (int i = 0; i < COLUMNS-1; i++) {
        for (int j = 0; j < ROWS; j++) {
            for (int k = 0; k < ROWS; k++) {
                if (i > 0) {
                    neurons[i][j]->output[k]->value += activationFunction(neurons[i][j]->value) * neurons[i][j]->weight[k];
                }
                else {
                    neurons[i][j]->output[k]->value += neurons[i][j]->value * neurons[i][j]->weight[k];
                }
            }
        }
    }
    
    return neurons[COLUMNS-1][0]->value;
}

int main () {
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Neural Network", sf::Style::Default, settings);
    sf::Image plotImage;

    // Create neurons
    for (int i = 0; i < COLUMNS; i++) {
        for (int j = 0; j < ROWS; j++) {
            neurons[i][j] = new Neuron;
            //printf("x ");
        }
        //printf("\n");
    }
    //printf("\n");

    // Intialise connections
    for (int i = 0; i < COLUMNS-1; i++) {
        for (int j = 0; j < ROWS; j++) {
            for (int k = 0; k < ROWS; k++) {
                neurons[i][j]->weight[k] = 0;
                neurons[i][j]->output[k] = neurons[i+1][k];
            }
        }
    }

    // Randomise population's genes
    for (int i = 0; i < POPULATION; i++) {
        for (int j = 0; j < COLUMNS-1; j++) {
            for (int k = 0; k < ROWS; k++) {
                for (int l = 0; l < ROWS; l++) {
                    genes[i][j][k][l] = (((float)rand() / (float)RAND_MAX)*2.0)-1.0;
                }
            }
        }
    }

    float last_loss = 0;
    int generation = 0;

    int best_members[POPULATION] = {0};
    float losses[POPULATION] = {INFINITY};

    // Run    
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // Clear screen
        plotImage.create(WIDTH/2, HEIGHT, sf::Color(0, 0, 0));
        window.clear(sf::Color(0,0,0,255));

        // Draw ideal output
        for (int x = 0; x < WIDTH/2; x++) {
            int y = (HEIGHT*0.5)+(targetFunction((((float)x / (float)WIDTH)*4.0)-1.0)*HEIGHT*0.5);
            plotImage.setPixel(x, y, sf::Color(255,255,255,255));
        }

        // Create new offspring to replace worst members
        for (int offspring = 0; offspring < OFFSPRING_PER_GENERATION; offspring++) {
            // Choose random points to splice genes between
            int crossover_low = (int)(((float)rand() / (float)RAND_MAX)*COLUMNS*ROWS*ROWS);
            int crossover_high = crossover_low + (int)(((float)rand() / (float)RAND_MAX)*((COLUMNS*ROWS*ROWS) - crossover_low));

            int parent_one = (int)(pow((float)rand() / (float)RAND_MAX, 2)*(float)(POPULATION-OFFSPRING_PER_GENERATION));
            int parent_two = (int)(pow((float)rand() / (float)RAND_MAX, 2)*(float)(POPULATION-OFFSPRING_PER_GENERATION));

            // Splice parent's genes to create offspring's genes
            for (int i = crossover_low; i < crossover_high; i++) {
                genes[best_members[POPULATION-offspring]][((i % ROWS)%ROWS)%COLUMNS][(i % ROWS)%ROWS][i % ROWS] = genes[best_members[parent_one]][((i % ROWS)%ROWS)%COLUMNS][(i % ROWS)%ROWS][i % ROWS];
            }
            for (int i = crossover_high; i >= crossover_low; i--) {
                genes[best_members[POPULATION-offspring]][((i % ROWS)%ROWS)%COLUMNS][(i % ROWS)%ROWS][i % ROWS] = genes[best_members[parent_two]][((i % ROWS)%ROWS)%COLUMNS][(i % ROWS)%ROWS][i % ROWS];
            }
        }

        // Add some tiny variations to the population
        for (int population = 0; population < POPULATION; population++) {
            for (int i = 0; i < COLUMNS-1; i++) {
                for (int j = 0; j < ROWS; j++) {
                    for (int k = 0; k < ROWS; k++) {
                        genes[best_members[population]][i][j][k] += pow((((float)rand() / (float)RAND_MAX)*(float)1.8)-(float)0.9, 25);
                    }
                }
            }
        }

        // Add a few larger mutations to the population biased towards worst in order to overcome local minima
        for (int mutation = 0; mutation < MUTATIONS; mutation++) {
            genes[best_members[POPULATION - (int)(pow((float)rand() / (float)RAND_MAX, 99)*(POPULATION))]][(int)(((float)rand() / (float)RAND_MAX)*(COLUMNS-1))][(int)(((float)rand() / (float)RAND_MAX)*ROWS)][(int)(((float)rand() / (float)RAND_MAX)*ROWS)] += pow((((float)rand() / (float)RAND_MAX)*(float)2.0)-(float)1.0, 15);
        }

        for (int population = 0; population < POPULATION; population++) {
            for (int i = 0; i < COLUMNS-1; i++) {
                for (int j = 0; j < ROWS; j++) {
                    for (int k = 0; k < ROWS; k++) {
                        genes[best_members[population]][i][j][k] = std::max(std::min(genes[population][i][j][k], (float)1), (float)-1);
                    }
                }
            }
        }

        // Simulate every member of population
        for (int population = 0; population < POPULATION; population++) {
            for (int i = 0; i < COLUMNS; i++) {
                for (int j = 0; j < ROWS; j++) {
                    for (int k = 0; k < ROWS; k++) {
                        neurons[i][j]->weight[k] = genes[population][i][j][k];
                    }
                }
            }
            // Run trials
            float loss = 0;
            for (int trial = 0; trial < TRIALS; trial++) {
                float input = (((float)trial / (float)TRIALS)*2.0)-1.0;
                float y = targetFunction(input);

                // Retrive output
                float output = simulate(input);

                // Punish big errors exponentially more than small ones
                loss += pow(fabs(y-output), 3);
            }

            losses[population] = loss;
            best_members[population] = population;
        }

        // Sort lists
        for (int i = 0; i < POPULATION-1; i++)      {   
            for (int j = 0; j < POPULATION-i-1; j++) {
                if (losses[j] < losses[j+1]) {
                    swap(&losses[j], &losses[j+1]);
                    swap(&best_members[j], &best_members[j+1]);
                }
            }
        }

        // Redraw screen for every generation
        if (true) {
            printf("%i - ", generation);
            printf("Best member: %i\t", best_members[0]);
            printf("Best loss: %f\t", losses[0]);
            printf("Improvement: %f\n", last_loss - losses[0]);

            // Set network to best of generation
            for (int i = 0; i < COLUMNS; i++) {
                for (int j = 0; j < ROWS; j++) {
                    for (int k = 0; k < ROWS; k++) {
                        neurons[i][j]->weight[k] = genes[best_members[0]][i][j][k];
                    }
                }
            }

            // Draw best output
            for (int x = 0; x < WIDTH/2; x++) {
                float input = (((float)x / (float)WIDTH)*4.0)-1.0;

                float output = simulate(input);
                output = std::min(std::max(output, (float)-1), (float)1);

                // Draw output
                int y = (HEIGHT*0.5)+(output*HEIGHT*0.5);
                plotImage.setPixel(x, y, sf::Color(255,0,0,255));
            }

            // Draw network
            for (int i = 0; i < COLUMNS; i++) {
                for (int j = 0; j < ROWS; j++) {
                    int x = (WIDTH/2)+((WIDTH/4)/COLUMNS)+(i*((WIDTH/2)/COLUMNS));
                    int y = (j*(HEIGHT/ROWS))+((HEIGHT/4)/ROWS);
                    int radius = HEIGHT/100;

                    if (i == 0) {
                        if (j > 0 and j < ROWS-1) {
                            continue;
                        }
                    }
                    if (i == COLUMNS-1) {
                        if (j > 0) {
                            continue;
                        }
                    }

                    sf::CircleShape circle (radius);
                    circle.setFillColor(sf::Color(255,255,255,255));
                    circle.setPosition(x, y);

                    for (int k = 0; k < ROWS; k++) {
                        if (i == COLUMNS-2) {
                            if (k > 0) {
                                continue;
                            }
                        }
                        int x1 = x+radius;
                        int y1 = y+radius;
                        int x2 = (WIDTH/2)+((WIDTH/4)/COLUMNS)+((i+1)*((WIDTH/2)/COLUMNS))+radius;
                        int y2 = (k*(HEIGHT/ROWS))+((HEIGHT/4)/ROWS)+radius;
                        int length = sqrt(pow((float)(x2-x1), 2) + pow((float)(y2-y1), 2));
                        float theta = (atan((float)(y2-y1)/(float)(x2-x1))/(2.0*M_PI))*360;

                        sf::RectangleShape line(sf::Vector2f(length, 1));
                        line.setPosition(x1, y1);
                        line.rotate(theta);
                        line.setFillColor(sf::Color((genes[best_members[0]][i][j][k] > 0.0)*255.0,
                                                    (genes[best_members[0]][i][j][k] < 0.0)*63.0,
                                                    (genes[best_members[0]][i][j][k] < 0.0)*255.0,
                                                    fabs(genes[best_members[0]][i][j][k]*127.0)));
                        window.draw(line);
                    }

                    window.draw(circle);
                }
            }

            sf::Texture texture;
            sf::Sprite sprite;
            texture.loadFromImage(plotImage);
            sprite.setTexture(texture, true);
            window.draw(sprite);
            window.display();

            last_loss = losses[0];
        }

        generation += 1;
    }
}